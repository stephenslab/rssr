function [A] = read_fvec(filename)
  fileId=fopen(filename,'r');
  formatspec='%f';
  A = fscanf(fileId,formatspec);
  fclose(fileId);
end
