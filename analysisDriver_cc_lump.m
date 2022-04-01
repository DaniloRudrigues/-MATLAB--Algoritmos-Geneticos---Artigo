function vpl = analysisDriver_cc_lump(dakotaResponseFile,var)

global ac

[tags] = ManipParamNp_cc_lump(var,ac);

templateFileName = ac.mxtpl;
runDir = ac.dir;
path(path,runDir)
templaterwd = ac.rrtpl;
imexFile = makeFileName_paralelo_vzs(runDir);

tempo = ac.time;

if tempo == 0
    processFilesNcte(templateFileName,fullfile(runDir, imexFile.Input), tags, ac,'\$');
else
    processFilesNvar_cc_lump(templateFileName,fullfile(runDir, imexFile.Input), tags, ac,'\$');
end

processRwd(templaterwd,fullfile(runDir,imexFile.rwd));       %novo!

cwd = pwd(); cd(runDir);

% imexCmd = ['mx201010.exe -f ' imexFile.Input ' -log'];     %
imexCmd = ['"C:\Program Files (x86)\CMG\IMEX\2015.10\Win_x64\EXE\mx201510.exe" -f ' imexFile.Input ' -log -wait'];
status = system(imexCmd);

if status ~= 0
  fprintf(1,'Error running imex ! \n');
%   exit(status);
end
reportCmd = ['"C:\Program Files (x86)\CMG\BR\2015.10\Win_x64\EXE\report.exe" -f ' imexFile.rwd '  -o ' imexFile.rwo];
status = system(reportCmd);                                         %
if status ~= 0                                                      % novo!
  fprintf(stderr,'Error running report ! \n');                      %
  exit(status);                                                     %
end                                                                 %

cd(cwd);

vpl = posProcessadorImexr(dakotaResponseFile,imexFile.rwo,ac);

cleanup(runDir, imexFile);

% printf('*****Finishing analysis driver ./n')