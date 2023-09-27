function writeResultsToTxt(fileOut, results, fileIn1, fileIn2)
% WRITERESULTSTOTXT Zapiši rezultate deformacijske analize z izbrano robustno metodo v datoteko.
%
% Inputs:
%    fileOut : ime izhodne datoteke
%    results : rezultati funkcije DARM
%    fileIn1 : ime vhodne datoteke prve  terminske izmere
%    fileIn2 : ime vhodne datoteke druge terminske izmere


    arguments
        fileOut string {mustBeTextScalar}
        results struct
        fileIn1 string {mustBeTextScalar}
        fileIn2 string {mustBeTextScalar}
    end
    
    wMethodsSloNames = configureDictionary('string', 'string');
    wMethodsSloNames("L1") = "metoda L1";
    wMethodsSloNames("L1-L2") = "metoda L1-L2";
    wMethodsSloNames("Lp") = "metoda Lp";
    wMethodsSloNames("Huber") = "Huberjeva metoda";
    wMethodsSloNames("HuberMod") = "modificirana Huberjeva metoda";
    wMethodsSloNames("Fair") = "Fairova metoda";
    wMethodsSloNames("Cauchy") = "Cauchyjeva metoda";
    wMethodsSloNames("Welsch") = "Welscheva metoda";
    wMethodsSloNames("Tukey") = "Tukeyeva metoda";
    wMethodsSloNames("Hampel") = "Hampelova metoda";
    wMethodsSloNames("German-McClure") = "German-McClurejeva metoda";
    wMethodsSloNames("danska") = "danska metoda";

    fout = fopen(fileOut, 'wt', 'n', 'UTF-8');
    
    fprintf(fout, '*** DEFORMACIJSKA ANALIZA Z ROBUSTNIMI METODAMI ***\n\n');
    fprintf(fout, 'datoteka prve terminske izmere : %s\n', fileIn1);
    fprintf(fout, 'datoteka druge terminske izmere: %s\n', fileIn2);
    fprintf(fout, '\n');
    fprintf(fout, 'Metoda dodeljevanja uteži: %s\n', wMethodsSloNames(results.wFunctionName));
    switch results.wFunctionName
        case 'Lp'
            fprintf(fout, '    * parameter c = %g', results.cParameter);
        case {'Huber', 'HuberMod', 'Fair', 'Cauchy', 'Welsch', 'Tukey', 'danska'}
            fprintf(fout, '    * parameter c = %.4f', results.cParameter);
        case 'Hampel'
            fprintf(fout, '    * parameter a = %.4f', results.aParameter);
            fprintf(fout, '    * parameter b = %.4f', results.bParameter);
            fprintf(fout, '    * parameter c = %.4f', results.cParameter);
    end
    fprintf(fout, '\n\n');
    
    logicalStr = ["            ", "      ✓     "];

    data = results.perComponents;
    fprintf(fout, '========================================================================================================================================\n');
    fprintf(fout, '================================================== ANALIZA STABILNOSTI PO KOMPONENTAH ==================================================\n');
    fprintf(fout, '========================================================================================================================================\n');
    fprintf(fout, '|     točka    |    dy [m]    |      Ty      | stabilnost y |    dx [m]    |      Tx      | stabilnost x |    ds [m]    |  stabilnost  |\n');
    fprintf(fout, '----------------------------------------------------------------------------------------------------------------------------------------\n');
    for i = 1 : numel(results.pointNames)
        fprintf(fout, '| %12s | %12.4f | %12.4f | %12s | %12.4f | %12.4f | %12s | %12.4f | %12s |\n', ...
            results.pointNames(i), data.d(2*i-1), data.tStatistic(i, 1), logicalStr(data.stableComponents(i, 1)+1), ...
                                   data.d(2*i  ), data.tStatistic(i, 2), logicalStr(data.stableComponents(i, 2)+1), ...
                                   data.s(i  ), logicalStr(data.stablePoints(i)+1));
    end
    fprintf(fout, '========================================================================================================================================\n');
    fprintf(fout, 'značilnost testa: %.2f %%\n', double(results.alpha)*100);
    fprintf(fout, 'T_krit          : %.4f\n', data.tCritValue);
    fprintf(fout, 'število iteracij: %d\n', data.nIter);
    fprintf(fout, '\n\n');

    
    data = results.perPoints;
    fprintf(fout, '===========================================================================================\n');
    fprintf(fout, '============================== ANALIZA STABILNOSTI PO TOČKAH ==============================\n');
    fprintf(fout, '===========================================================================================\n');
    fprintf(fout, '|     točka    |    dy [m]    |    dx [m]    |    ds [m]    |      T       |  stabilnost  |\n');
    fprintf(fout, '-------------------------------------------------------------------------------------------\n');
    for i = 1 : numel(results.pointNames)
        fprintf(fout, '| %12s | %12.4f | %12.4f | %12.4f | %12.4f | %12s |\n', ...
            results.pointNames(i), data.d(2*i-1), data.d(2*i), data.s(i), data.tStatistic(i), logicalStr(data.stablePoints(i)+1));
    end
    fprintf(fout, '===========================================================================================\n');
    fprintf(fout, 'značilnost testa: %.2f %%\n', double(results.alpha)*100);
    fprintf(fout, 'T_krit          : %.4f\n', data.tCritValue);
    fprintf(fout, 'število iteracij: %d\n', data.nIter);
    fprintf(fout, '\n\n');
    
    fclose(fout);
end

