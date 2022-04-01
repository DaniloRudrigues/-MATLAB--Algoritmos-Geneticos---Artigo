function processFilesNvar_cc_lump(templateFile, outputFile, tags, ac, marker)
%
%  processFiles: 
%    Copy template file to output file, substituting the
%    variables indicated in the tags array by their
%    corresponding values.
%    templateFile: (string) name of input template file.
%    outputFile: (string) name of processed output file.
%    tags: (array of structure) name/value pairs.
%    marker: (string) delimiter for variable

%   This will scan the file looking for variables, register
%   where each variable begins, then write out to the output
%   file a segment where there are no variables, write the value of
%   the variable, write the next segment, and so on until the end of
%   file.

    n_interv_tempo = ac.ncc;
    Tc = ac.T; % Tempo de concessão do reservatório
    ano = 2006;		    % Data inicial 2006 01 01

    [inpf, msg] = fopen(templateFile, 'rb');
    if inpf == -1
        fprintf(1, 'Error opening file %s\n', templateFile);
        fprintf(1, 'Error message is: ', msg);
%         exit(1)
    end

    [outf, msg] = fopen(outputFile, 'wb');
    if outf == -1
        fprintf(1, 'Error opening file %s\n', outputFile);
        fprintf(1, 'Error message is: ', msg);
%         exit(1)
    end

    [fileContents, count] = fscanf(inpf, '%c', Inf);
    
    % Create a pattern to define a variable.
    % Starts with the marker, followed by one of more word characters,
    % and ends with a non word character.
    varPattern = [ marker, '\w+' ];

    % Find out where all variables start and end.
    [vstart, vend, vname ] = ...
        regexp(fileContents, varPattern, 'start', 'end', 'match');

    for i = 1: length(tags) % O octave não aceitou classificar tipo struct
        aux = tags(i).time;
        ch_time(i,1) = round(aux);
        ch_time(i,2) = i;
    end
    
    % Cria um vetor de tempos
   
    time = calendario(ac,ano);
    [lt, ct] = size(time);

    % Inclui no vetor de tempo onde ocorrerá mudança

    for i = 1:size(ch_time,1)
        int = lt + i;
        time(int,1) = ch_time(i,1);
        time(int,2) = ch_time(i,2);
    end
    time = sortrows(time,1);

    % Remove intervalos, não marcados para alterações, do vetor que tenha intervalos menores que um dia	
    
    remo = [];
    y = diff(time(:,1));
    menorq = find(y<1);
    ii = 1;

    for i=1:length(menorq)
        if (time(menorq(i),2)==0)
            remo(ii) = menorq(i);
            ii = ii + 1;
        elseif (time(menorq(i)+1,2)==0)
            remo(ii)= menorq(i)+1;
            ii = ii + 1;
        end
    end
    if length(remo)>0
        time(remo,:) = [];
    end
	
    % Write output file, replacing times, rates for wells and variables by their values .
	
    n_pocos_prod = ac.npp;
    n_pocos_inje = ac.npi;
    segStart = 1;

    fprintf(outf, '%s',fileContents(segStart:vstart(1)-1));
    
    for it = 1:length(time)

	if it > 1  
		s_impre = time(it,1)-time(it-1,1);	
		if s_impre > 0		
        	impre = num2str(time(it,1));
            fprintf(outf, '%s', '*TIME ', impre);
            fprintf(outf, '\n\n');
		end
	end
%%%%%%%%%%%%%%%%%%%%%%%% MEXI AQUI %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%	
	if time(it,2)>0
		fprintf(outf, '%s', ' *CLUMPSETTING ');
        fprintf(outf, '\n');
		if tags(time(it,2)).type == 1
            impre = sprintf('''LUMP_PROD%i''', tags(time(it,2)).number);
            fprintf(outf, '%s', impre);
% 		elseif tags(time(it,2)).type == 2
%             impre = sprintf('''INJ%i''', tags(time(it,2)).number);
%             fprintf(outf, '%s', impre);
		end
            fprintf(outf, '\t');				     		
        	fprintf(outf, '%g', tags(time(it,2)).val);
            fprintf(outf, '\t');
            fprintf(outf, '\n\n');
	end
    end	
    % Last segment
    segStart = vend(1) + 1;
    % Last segment
    fprintf(outf, '%s',fileContents(segStart:end));

    fclose(inpf);
    fclose(outf);
