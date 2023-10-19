  function prepend2file( string, filename, newline )
% function prepend2file( string, filename, newline )
%
% 
%  newline: is an optional boolean, that if true will append a \n to the end
%  of the string that is sent in such that the original text starts on the
%  next line rather than the end of the line that the string is on
%  string: a single line string
%  filename: the file you want to prepend to
%
% The function prepend2file is part of the MATLAB toolbox SciXMiner. 
% Copyright (C) 2010  [Ralf Mikut, Tobias Loose, Ole Burmeister, Sebastian Braun, Andreas Bartschat, Johannes Stegmaier, Markus Reischl]


% Last file change: 15-Aug-2016 11:31:14
% 
% This program is free software; you can redistribute it and/or modify,
% it under the terms of the GNU General Public License as published by 
% the Free Software Foundation; either version 2 of the License, or any later version.
% 
% This program is distributed in the hope that it will be useful, but
% WITHOUT ANY WARRANTY; without even the implied warranty of 
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License along with this program;
% if not, write to the Free Software Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110, USA.
% 
% You will find further information about SciXMiner in the manual or in the following conference paper:
% 
% 
% Please refer to this paper, if you use SciXMiner for your scientific work.

tempFile = tempname
fw = fopen( tempFile, 'wb' );
if nargin < 3
newline = false;
end
if newline
fwrite( fw, sprintf('%s\n', string ) );
else
fwrite( fw, string );
end

fclose( fw );
appendFiles( filename, tempFile );
copyfile( tempFile, filename );
delete(tempFile);

% append readFile to writtenFile
function status = appendFiles( readFile, writtenFile )
fr = fopen( readFile, 'rb' );
fw = fopen( writtenFile, 'ab' );

while feof( fr ) == 0
tline = fgetl( fr );
fwrite( fw, sprintf('%s\n',tline ) );
end
fclose(fr);
fclose(fw);