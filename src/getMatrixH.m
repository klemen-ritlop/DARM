function H = getMatrixH(yxPoints, observedDistances)
% GETMATRIXH Sestavi matriko H.
%
% Inputs:
%    yxPoints          : matrika s koordinatami točk (vsaka točka nova vrstica, prvi
%                        stolpec koordinate y, drugi stolpec koordinate x)
%    observedDistances : ali imamo v geodetski mreži opazovane dolžine (*true* / false)
%
% Outputs:
%    H                 : matrika H

    arguments
        yxPoints (:,2) {mustBeNumeric, mustBeReal}
        observedDistances logical {mustBeNumericOrLogical} = true
    end
    
    nPoints = size(yxPoints, 1);
    y = yxPoints(:, 1);
    x = yxPoints(:, 2);

    % koordinate reduciramo na težišče mreže
    y = y - mean(y);
    x = x - mean(x);
    
    % s H_ je označeno H'
    H_ = zeros(4, 2*nPoints);
    H_(1, 1:2:end-1) = ones(1, nPoints);
    H_(2, 2:2:end)   = ones(1, nPoints);
    H_(3, 1:2:end-1) = x;
    H_(3, 2:2:end)   = -y;
    H_(4, 1:2:end-1) = y;    
    H_(4, 2:2:end)   = x;
    
    % vsako vrstico matrike H delimo z normo pripadajoče vrstice (numerično bolj stabilno)
    for i = 1 : size(H_, 1)
        H_(i, :) = H_(i, :) / norm(H_(i, :));
    end
    
    % če imamo merjene dolžine, izpustimo zadnjo vrstico
    if observedDistances
        H_ = H_(1:3, :);
    end
    
    % s H_ je označeno H'
    H = H_';
    
end

