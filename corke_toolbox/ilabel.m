%ILABEL	Label an image
%
%		li = ilabel(image [, connect])
%		[li,maxlabel] = ilabel(image [, connect])
%		[li,maxlabel,parents] = ilabel(image [, connect])
%
%	Perform connectivity analysis and return a label image LI, same size as 
%	IMAGE where each pixel value represents the label of the region of its
%	corresponding pixel in IMAGE.
%
%	Connectivity is performed using 4 or 8 nearest neighbour as controlled
%	by the optional connect argument.  Default is 4 way but can be changed
%	to 8.  The function compares the pixel values so the image can be 
%	greyscale or binary.
%
%	Also returns the maximum label assigned to the image.  Labels lie in
%	the range 1 to MAXLABEL.
%
%	The third form also returns region hierarchy information.  The
%	value of parents[i] is the label of the 'parent' or enclosing
%	region of region i.  A value of 0 indicates that the region has
%	no enclosing region.
%
% SEE ALSO:	imoments iblobs
%
% LIMITATIONS:	should allow for different connectivity modes
%
%	  By Peter Corke, Copyright (c) CSIRO, 1999  Machine Vision Toolbox for Matlab

% $Header: /home/autom/pic/cvsroot/image-toolbox/ilabel.m,v 1.4 2005/11/24 11:14:36 pic Exp $
% $Log: ilabel.m,v $
% Revision 1.4  2005/11/24 11:14:36  pic
% Change to copyright notices.
%
% Revision 1.3  2005/10/30 03:28:00  pic
% Doco tidyup
%
% Revision 1.2  2005/07/03 10:47:52  pic
% Add support for 4 and 8-way connectivity.
%
% Revision 1.1.1.1  2002/05/26 10:50:22  pic
% initial import
%

