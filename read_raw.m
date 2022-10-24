function vid = read_raw(filename, pix_w, pix_h, varargin)
% FUNCTION NAME:
%   read_raw
%
% DESCRIPTION:
%   Reads video file into an array
%trash
% INPUT (REQUIRED):
%        filename: (string) filename of video to be imported
%           pix_w: (double) pixel width of video frames
%           pix_h: (double) pixel height of video frames
%
% INPUT (OPTIONAL):
%     byte_offset: (double) Offset to first image (bytes to skip)
%    byte_spacing: (double) Gap between images (bytes to skip)
%     frame_start: (double) Starting frame to import (default is 1)
%         frame_N: (double) Number of frames to import (detault is all)
%      frame_skip: (double) Frames to skip (1 will skip every other frame)
%          packed: set to 'y' if data is 10-bit packed.
%        unpacked: set to 'y' if data is 10-bit unpacked.
%          max255: set to 'y' to output video to 8-bit
%
% OUTPUT:
%             vid: Output video (pix_h, pix_w, frame_N)
%
% CALLING SEQUENCE:
%   a = read_raw(filename, pix_w, pix_h)
%   a = read_raw('test.vif',656,491,byte_offset=72,byte_spacing=136,unpacked='y')
%   a = read_raw('test.vif',700,700,byte_offset=8,byte_spacing=496, frame_start=100, frame_N=200)
%
% NOTES :
%   ImageJ - ImageJ (free software) is useful for determining the
%           parameters needed to read in a raw video format. With ImageJ
%           open, use File->Import->Raw.  Playaround with parameters as
%           needed.
%   VIF DataType - XCAP software saves video in a VIF format using three
%           possible data types: 8-bit unsigned, packed (10-bit unsigned
%           packed into bytes, where 5 bytes [40 bits] is 4 data points [40
%           bits]), and unpacked (10-bit unsigned interger saved in 16-bit
%           unsigned interger, so 6 bits are unused).
%   VIF pix_w,pix_h - In the .fmt file, saved with the VIF file,
%           is information about these values (xviddim and yviddim).
%   VIF byte_spacing - To determine byte_spacing, use the the frameIOSize
%           in the .ini file. To determine byte_spacing, subtract number of
%           bytes for each image from this vale. Example: 700x700 8-bit (1
%           byte) image with a frameIOsize of 490496 will result in
%           byte_spacing = 490496-700*700*1 byte or 496 bytes.
%
% REVISION HISTORY:
%   05/30/2013 - K Aptowicz
%       * Wrote original version in IDL
%   10/23/2022 - K Aptowicz
%       * Translated to MATLAB
%
%% Reading and setting parameters
% Set default values for optional parameters
default_byte_offset = 0;
default_byte_spacing = 0;
default_frame_start = 1;
default_frame_skip = 0;
default_frame_N = [];
default_packed = [];
default_unpacked = [];
default_max255 = [];

% Create fields for all optionals inputs
%Variables
p = inputParser;
addParameter(p,'byte_offset',default_byte_offset,@isnumeric)
addParameter(p,'byte_spacing',default_byte_spacing,@isnumeric)
addParameter(p,'frame_start',default_frame_start,@isnumeric)
addParameter(p,'frame_skip',default_frame_skip,@isnumeric)
addParameter(p,'frame_N',default_frame_N,@isnumeric)
% Keywords
addOptional(p,'packed',default_packed)
addOptional(p,'unpacked',default_unpacked)
addOptional(p,'max255',default_max255)

% populate optional parameters from inputs
parse(p,varargin{:})
byte_offset = p.Results.byte_offset;
byte_spacing = p.Results.byte_spacing;
frame_start = p.Results.frame_start;
frame_skip = p.Results.frame_skip;
frame_N = p.Results.frame_N;
packed = p.Results.packed;
unpacked = p.Results.unpacked;
max255 = p.Results.max255;

% Image size
m = pix_w;
n = pix_h;

% Read in file
finfo=dir(filename);
if size(finfo,1) == 0
    disp(['No files matched specification ',filespec])
else
    fid = fopen(finfo(1).name);
end

% Determine maximum number of frames
Nbytes=finfo.bytes;
if ~isempty(packed)
    disp('Reading a packed video file was never tested!')
    if mod(m*n*5/4,1) ~= 0
        disp('Issue with bitpacking. pix_w*pix_h*10bit/8bit result in non-interger number of bytes.')
    end
    fmax=floor((Nbytes-byte_offset-m*n*5/4)./(m*n*5/4+byte_spacing))+1;
    fsize=m*n*5/4+byte_spacing
    data_type = 'uint16';
elseif ~isempty(unpacked)
    fmax=floor((Nbytes-byte_offset-m*n*2)./(m*n*2+byte_spacing))+1;
    data_type = 'uint16';
    fsize=m*n*2+byte_spacing;
else
    fmax=floor((Nbytes-byte_offset-m*n)./(m*n+byte_spacing))+1;
    data_type = 'uint8';
    fsize=m*n+byte_spacing;
end

% Determine last frame number (fmax) and number of frames
if ~isempty(frame_N)
    frame_end = frame_start + frame_N;
    frame_N = floor((frame_N-1)/(frame_skip+1))+1;
else
    frame_N = fmax-frame_start+1;
    frame_end=fmax;
    frame_N = floor((frame_N-1)/(frame_skip+1))+1;
end

% Goto first frame to be read in
fseek(fid, byte_offset+int32(frame_start-1)*fsize, 'bof');

%% Read in frames
vid=zeros(n,m,frame_N,data_type); % pre-allocation for video
for count=1:frame_N
    % Skip byte spacing between frames
    if ~isempty(byte_spacing) && count > 1
        fseek(fid,byte_spacing,'cof'); % Remove extra bytes between frames
    end
    % Skip frames if requested
    if (frame_skip~=0) && count > 1
        fseek(fid,fsize*frame_skip,'cof'); % Remove extra bytes between frames
    end
    if ~isempty(packed)
        frame = fread(fid,int64(m*n*5/4),'ubit10->uint16');

    else
        rawimg=fread(fid,[m,n],data_type);
    end
    vid(:,:,count)=rawimg';
    count=count+1;
end

if ~isempty(max255) && (data_type ~= "uint8")
    vid=uint8(vid/4);
end

fclose(fid);

end

%
% ;+
% ; NAME:
% ;      READ_RAW
% ;
% ; PURPOSE:
% ;      Read in video files which are a sequence of images.
% ;
% ; CATEGORY:
% ;      Kevin/Utility
% ;
% ; CALLING SEQUENCE:
% ;      a = read_raw(filename, pix_w, pix_h)
% ;      a = read_raw('test.vif',656,491,byte_offset=72,byte_spacing=136,unpacked='y')
% ;      a = read_raw('test.vif',700,700,byte_offset=8,byte_spacing=496, frame_start=100, frame_N=200)
% ;
% ; INPUTS:
% ;       filename: name of binary video file to be opened
% ;          pix_w: width of frame in pixels
% ;          pix_h: width of frame in pixels .
% ;    byte_offset: number of bytes to skip before first image
% ;   byte_spacing: number of bytes to skip between frames
% ;
% ; KEYWORDS:
% ;    frame_start: number of frames to skip before reading 'first' frame
% ;       frame_N: total number of frames to read in
% ;     frame_skip: Number of frames to skip (e.g. 1 will skip every other frame)
% ;         packed: set this keyword if data is 10 bit packed.
% ;       unpacked: set this keyword if data is 10 bit unpacked.
% ;          scale: set this keyword to scale the outmatrix into 8 bits
% ;
% ; OUTPUTS:
% ;              a: three dimensional matrix containing video
% ;
% ; RESTRICTIONS:
% ;
% ; PROCEDURE:
% ;      Wrote to be similar to ImageJ import function.
% ;
% ; NOTES:
% ;    KBA: To use to read in .vif files from uniqvision packed 10-bit CCD, the following worked
% ;         IDL> a=read_raw('frames_10_1fps.vif',656,491,byte_offset=72,byte_spacing=324,/packed, frame_N=10)
% ;
% ;    KBA: To determine input parameters, use the .ini file
% ;         frameIOSize --> Bytes for each frame. To determine byte_spacing, subtract number of bytes for each image
% ;         example ... a 700 x 700 8-bit (1 byte) image with a frameIOsize of 490496 will result in
% ;         byte_spacing = 490496-700*700*1 byte
% ;
% ;
% ; MODIFICATION HISTORY:
% ; Created by Kevin B. Aptowicz, West Chester University
% ; 30-May-2013: KBA.
%
% ;
% ; Copyright (c) 2013 Kevin B. Aptowicz
% ;
% ;-

% ; Get the house in order ...
% m = long64(pix_w)
% n = long64(pix_h)
% if  not keyword_set(byte_offset) then byte_offset=0
% if  not keyword_set(byte_spacing) then byte_spacing=0
% if  not keyword_set(frame_start) then frame_start=0
%
% ; Determine maximum number of frames
% finfo=file_info(filename)
% Nbytes=finfo.size
% if keyword_set(packed) then begin
% 	fmax=floor((Nbytes-byte_offset-m*n*5L/4L)/(m*n*5L/4L+byte_spacing))+1
% 	rawimg = bytarr(n*m*5L/4L)              ;; 10-bit packed is 20% more than 8-bit
% 	fsize=m*n*5L/4L+byte_spacing
% endif else if keyword_set(unpacked) then begin
%     fmax=floor((Nbytes-byte_offset-m*n)/(m*n*2L+byte_spacing))+1
%     rawimg = uintarr(n*m)              ;; 10-bit unpacked is 100% more than 8-bit
%     fsize=m*n*2L+byte_spacing
% endif else begin
%     fmax=floor((Nbytes-byte_offset-m*n)/(m*n+byte_spacing))+1
%     rawimg = bytarr(n*m)              ;; plane old 8-bit
%     fsize=m*n+byte_spacing
% endelse
%
% if keyword_set(frame_N) then begin
% 	if fmax lt (frame_start + frame_N) then begin
% 	print, 'WARNING: Not enough frames ... truncating size'
% 	if frame_start gt fmax   then frame_start=0
% 	endif else begin
% 		fmax = frame_start + frame_N
% 	endelse
%
% endif
%
%
%
% vid=uintarr(m,n,fmax-frame_start)
% count=0L
% for i=frame_start, fmax-1 do begin
% 	if keyword_set(byte_spacing) && count ne 0 then readu,lun,trash
% 	readu,lun, rawimg             ;; Read data
%     if keyword_set(packed) then begin
%     	rawimg = reform(rawimg, 5L, n*m/4L)     ;; Reform into quintuples
% 		rawimg=reverse(rawimg)
% 		outimg = intarr(4, n*m/4)             ;; Make new output image of quadruples
% 		outimg(0,*) = ishft(rawimg(0,*) AND 'ff'x,2) + ishft(rawimg(1,*) AND 'c0'x,-6)
% 		outimg(1,*) = ishft(rawimg(1,*) AND '3f'x,4) + ishft(rawimg(2,*) AND 'f0'x,-4)
% 		outimg(2,*) = ishft(rawimg(2,*) AND '0f'x,6) + ishft(rawimg(3,*) AND 'fc'x,-2)
% 		outimg(3,*) = ishft(rawimg(3,*) AND '03'x,8) + ishft(rawimg(4,*) AND 'ff'x, 0)
% 		outimg=reverse(outimg)
% 		outimg = reform(outimg, n, m)         ;; Convert back to an NxM array
% 		vid(*,*,count)=outimg
% 	endif else if keyword_set(unpacked) then begin
% 	vid(*,*,count)=reform(rawimg,n,m)
% 	endif else vid(*,*,count)=reform(rawimg,n,m)
% 	count=count+1L
% endfor
%
% if keyword_set(scale) then begin
% 	vid=double(vid)
% 	vid=byte(vid/4L)
% endif
%
% free_lun, lun
% return, vid
% end