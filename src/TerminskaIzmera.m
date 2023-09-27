classdef TerminskaIzmera
    % RAZRED TERMINSKAIZMERA
    %   V razred TerminskaIzmera shranimo podatke, ki so potrebni
    %   za izvedbo deformacijske analize z robustnimi metodami.
    %   Trenutno je implementirano branje podatkov iz datoteke *.daf,
    %   možna pa je enostavna implementacija branja iz drugih formatov.
    
    properties
        points
        nPoints
        redundancyNumber
        datumDefect
        s0Apost
        matrixQ
        matrixN
    end
    
    methods
        function obj = TerminskaIzmera( ...
                pointNames, ...
                yApproxCrd, ...
                xApproxCrd, ...
                yAdjCrd, ...
                xAdjCrd, ...
                redundancyNumber, ...
                datumDefect, ...
                s0Apost, ...
                matrixQ, ...
                matrixN ...
                )
            
            arguments
                pointNames {mustBeVector, mustBeText}
                yApproxCrd {mustBeVector, mustBeReal}
                xApproxCrd {mustBeVector, mustBeReal}
                yAdjCrd {mustBeVector, mustBeReal}
                xAdjCrd {mustBeVector, mustBeReal}
                redundancyNumber {mustBeScalarOrEmpty, mustBeNonempty, mustBeReal, mustBePositive}
                datumDefect {mustBeScalarOrEmpty, mustBeNonempty, mustBeReal, mustBeNonnegative}
                s0Apost {mustBeScalarOrEmpty, mustBeNonempty, mustBeReal, mustBeNonnegative}
                matrixQ {mustBeNumeric, mustBeReal}
                matrixN {mustBeNumeric, mustBeReal}
            end


            nPoints = length(pointNames);
            

            if length(yApproxCrd) ~= nPoints
                error('%s: Dimenzija vektorja yApproxCrd ni skladna z dimenzijo vektorja pointNames.', mfilename("class"))
            end

            if length(xApproxCrd) ~= nPoints
                error('%s: Dimenzija vektorja xApproxCrd ni skladna z dimenzijo vektorja pointNames.', mfilename("class"))
            end

            if length(yAdjCrd) ~= nPoints
                error('%s: Dimenzija vektorja yAdjCrd ni skladna z dimenzijo vektorja pointNames.', mfilename("class"))
            end
            
            if length(xAdjCrd) ~= nPoints
                error('%s: Dimenzija vektorja xAdjCrd ni skladna z dimenzijo vektorja pointNames.', mfilename("class"))
            end

            if size(matrixQ, 1) ~= size(matrixQ, 2)
                error('%s: Matrika Q mora biti kvadratna matrika.', mfilename("class"))
            end

            if size(matrixQ, 1) ~= 2*nPoints || size(matrixQ, 2) ~= 2*nPoints
                error('%s: Dimenzije matrike Q ne ustrezajo številu točk.', mfilename("class"))
            end

            if size(matrixN, 1) ~= size(matrixN, 2)
                error('%s: Matrika N mora biti kvadratna matrika.', mfilename("class"))
            end

            if size(matrixN, 1) ~= 2*nPoints || size(matrixN, 2) ~= 2*nPoints
                error('%s: Dimenzije matrike N ne ustrezajo številu točk.', mfilename("class"))
            end
            
            
            obj.nPoints = length(pointNames);
            obj.redundancyNumber = redundancyNumber;
            obj.datumDefect = datumDefect;
            obj.s0Apost = s0Apost;
            obj.matrixQ = matrixQ;
            obj.matrixN = matrixN;

            obj.points = struct( ...
                'name', cellstr(pointNames), ...
                'yApprox', num2cell(yApproxCrd), ...
                'xApprox', num2cell(xApproxCrd), ...
                'yAdj', num2cell(yAdjCrd), ...
                'xAdj', num2cell(xAdjCrd) ...
                );
            
            % pretvori imena ('name') iz char array v string
            tmp = num2cell(string({obj.points.name}));
            [obj.points.name] = tmp{:};

        end
    end
        
        
    methods(Static)
        % --------------------------------------------------------------
        % KONSTRUKTOR, KI NAREDI OBJEKT NA PODLAGI VHODNE *.daf DATOTEKE
        % --------------------------------------------------------------
        function obj = fromDafFile(dafFile)
                arguments
                    dafFile {mustBeFile}
                end
                
                [~, name, ext] = fileparts(dafFile);
                dafFileName = [name, ext];

                % zastavice v vhodni datoteki
                fl  = '*';     % kontrolni znak
                flK = '*K';    % konec datoteke
                flF = '*F';    % število nadštevilnih merjenj (po izravnavi)
                flT = '*T';    % izravnane koordinate točk
                flR = '*R';    % defekt datuma
                flM = '*M';    % a-posteriori srednji pogrešek utežne enote
                flQ = '*Q';    % matrika kofaktorjev koordinatnih razlik
                flN = '*N';    % matrika normalnih enačb
                
            
                finp = fopen(dafFile, 'rt');
            
                % pridobi število točk v vhodni datoteki
                nPoints_ = 0;
                line = fgetl(finp);
                stop = false;
                while ischar(line) && ~stop
                    if length(line) < 2
                        line = fgetl(finp);
                        continue
                    end
                    
                    if strcmp(line(1:2), flT)
                        while ischar(line)
                            line = fgetl(finp);
                            if length(line) < 2
                                continue
                            elseif strcmp(line(1), fl)
                                stop = true;
                                break
                            else
                                nPoints_ = nPoints_ + 1;
                            end
                        end
                    else
                        line = fgetl(finp);
                    end
                    
                end
            
                if nPoints_ == 0
                    error('Razred %s: V datoteki %s ni točk.', mfilename("class"), dafFileName)
                end
            
            
                % nastavi dimenzije spremenljivk
                pointNames = strings(nPoints_, 1);        % imena točk
                approxCrd = zeros(nPoints_, 2);           % približne koordinate točk
                adjCrd = zeros(nPoints_, 2);              % izravnane koordinate točk
                Qv = zeros((2*nPoints_)^2, 1);            % matrika kofaktorjev koordinatnih razlik (vse vrstice zložene v en vektor)
                Nv = zeros((2*nPoints_)^2, 1);            % matrika normalnih enačb (vse vrstice zložene v en vektor)
                
            
                % preberi datoteko s podatki
                frewind(finp)
                line = fgetl(finp);
                rowCounterN = 0;
                rowCounterQ = 0;
                rowCounterT = 0;
            
                while ischar(line)                % preskoči prazne vrstice
                    if isempty(line)
                        line = fgetl(finp);
                        continue
                    end
                    
                    if strcmp(line(1), fl)        % zastavica za konec datoteke -> končaj branje
                        if strcmp(line(1:2), flK)
                            break
                        else
                            section = line(1:2);
                            line = fgetl(finp);
                            continue
                        end
                    end
            
                    switch section
                        % preberi število nadštevilnih merjenj (po izravnavi)
                        case flF
                            redundancyNumber_ = str2double(line);
            
                        % preberi koordinate točk
                        case flT
                            line = replace(line, "'", " ");
                            line = strsplit(strip(line));
                            pointNames(rowCounterT + 1) = string(line{1});
                            approxCrd(rowCounterT + 1, :) = [str2double(line{2}), str2double(line{3})];
                            adjCrd(rowCounterT + 1, :) = [str2double(line{4}), str2double(line{5})];
                            rowCounterT = rowCounterT + 1;
                        
                        % preberi defekt datuma
                        case flR
                            datumDefect_ = str2double(line);
                        
                        % preberi a-posteriori srednji pogrešek utežne enote
                        case flM
                            s0Apost_ = str2double(line);
            
                        % preberi matriko kofaktorjev koordinatnih razlik
                        case flQ
                            line = str2double(strsplit(strip(line)))';
                            nelements = length(line);
                            rowCounterQ = rowCounterQ + nelements;
                            Qv(rowCounterQ-nelements+1:rowCounterQ, 1) = line;
            
                        % preberi matriko normalnih enačb
                        case flN
                            line = str2double(strsplit(strip(line)))';
                            nelements = length(line);
                            rowCounterN = rowCounterN + nelements;
                            Nv(rowCounterN-nelements+1:rowCounterN, 1) = line;
                    end
                    
                    line = fgetl(finp);
            
                end
            
                fclose(finp);
                
                % preveri prebrane podatke
                if sum(Qv) == 0
                    error('Razred %s: Datoteka %s - matrika Q ni bila uspešno prebrana.', mfilename("class"), dafFileName)
                end
            
                if sum(Nv) == 0
                    error('Razred %s: Datoteka %s - matrika N ni bila uspešno prebrana.', mfilename("class"), dafFileName)
                end
            
                if rowCounterQ ~= (2*nPoints_)^2
                    error('Razred %s: Datoteka %s - prebrana matrika Q ni ustreznih dimenzij.', mfilename("class"), dafFileName)
                end
            
                if rowCounterN ~= (2*nPoints_)^2
                    error('Razred %s: Datoteka %s - prebrana matrika N ni ustreznih dimenzij.', mfilename("class"), dafFileName)
                end

                % preoblikuj vektor Qv v kvadratno matriko Q
                matrixQ_ = reshape(Qv, 2*nPoints_, 2*nPoints_)';
            
                % preoblikuj vektor Nv v kvadratno matriko N
                matrixN_ = reshape(Nv, 2*nPoints_, 2*nPoints_)';

                obj = TerminskaIzmera( ...
                    pointNames, ...
                    approxCrd(:,1), ...
                    approxCrd(:,2), ...
                    adjCrd(:,1), ...
                    adjCrd(:,2), ...
                    redundancyNumber_, ...
                    datumDefect_, ...
                    s0Apost_, ...
                    matrixQ_, ...
                    matrixN_ ...
                    );
        end
    end
end

