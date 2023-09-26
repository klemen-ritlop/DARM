function W = getMatrixW(d, dSigma, wFunctionName, a, b, c)
% GETMATRIXW Sestavi matriko uteži W za S-transformacijo.
%
% Inputs:
%    d             : vektor premikov
%    dSigma        : vektror standardnih odklonov prmikov v vektorju d
%    wFunctionName : ime utežne funkcije (*L1*, L1-L2, Lp, Huber, HuberMod, Fair,
%                    Cauchy, Welsch, Tukey, German-McClure, Hampel, danska)
%    a, b, c       : parametri nekaterih utežnih funkcij
%
% Outputs:
%     W            : matrika uteži W za S-transformacijo


    arguments
        d {mustBeVector, mustBeNumeric, mustBeReal}
        dSigma {mustBeVector, mustBeNumeric, mustBeReal}
        wFunctionName string {mustBeMember(wFunctionName, {'L1', 'L1-L2', 'Lp', 'Huber', 'HuberMod', ...
            'Fair', 'Cauchy', 'Welsch', 'Tukey', 'German-McClure', 'Hampel', 'danska'})} = 'L1'
        a {mustBeNumeric, mustBeReal} = -1
        b {mustBeNumeric, mustBeReal} = -1
        c {mustBeNumeric, mustBeReal} = -1
    end
    

    % privzete vrednosti parametrov posameznih utežnih funkcij
    switch wFunctionName
        case 'Lp'
            if c == -1; c = 1.2000; end
        case 'Huber'
            if c == -1; c = 1.3450; end
        case 'HuberMod'
            if c == -1; c = 1.2107; end
        case 'Fair'
            if c == -1; c = 1.3998; end
        case 'Cauchy'
            if c == -1; c = 2.3849; end
        case 'Welsch'
            if c == -1; c = 2.9846; end
        case 'Tukey'
            if c == -1; c = 4.6851; end
        case 'danska'
            if c == -1; c = 3.0000; end
        case 'Hampel'
            if a == -1 && b == -1 && c == -1
                a = 1.5000;
                b = 3.0000;
                c = 6.0000;
            end
    end
    

    n = length(d);
    W_ = zeros(n, 1);    % W_ - vektor diagonalnih elementov matrike W
    
    switch wFunctionName
        
        case 'L1'
            W_ = 1./abs(d);

        case 'L1-L2'
            W_ = 1 ./ sqrt(1 + d.^2 ./ 2);

        case 'Lp'
            nu = c;
            W_ = abs(d).^(nu-2);

        case 'Huber'
            q = c * dSigma;
            cond = abs(d) <= q;
            W_(cond) = 1;
            W_(~cond) = q(~cond) ./ abs(d(~cond));

        case 'HuberMod'
            q = c * dSigma;
            cond = (abs(d) ./ q) <= pi/2;
            W_(cond) = q(cond) ./ abs(d(cond) + 1e-16) .* sin(abs(d(cond) + 1e-16) ./ q(cond));
            W_(~cond) = q(~cond) ./ abs(d(~cond));

        case 'Fair'
            q = c * dSigma;
            W_ = 1 ./ (1 + abs(d) ./ q);
        
        case 'Cauchy'
            q = c * dSigma;
            W_ = 1 ./ (1 + (abs(d) ./ q).^2);
        
        case 'Welsch'
            q = c * dSigma;
            W_ = exp(-(abs(d) ./ q).^2);

        case 'Tukey'
            q = c * dSigma;
            cond = abs(d) <= q;
            W_(cond) = (1 - (d(cond) ./ q(cond)).^2).^2;
            W_(~cond) = 0;

        case 'German-McClure'
            W_ = 1 ./ (1+(d.^2)).^2;

        case 'Hampel'
            q = a * dSigma;
            u = b * dSigma;
            v = c * dSigma;
            cond1 = abs(d) <= q;
            cond2 = abs(d) > q & abs(d) <= u;
            cond3 = abs(d) > u & abs(d) <= v;
            cond4 = abs(d) > v;
            W_(cond1) = 1;
            W_(cond2) = q(cond2) ./ abs(d(cond2));
            W_(cond3) = q(cond3) .* (v(cond3) - abs(d(cond3))) ./ (abs(d(cond3)) .* (v(cond3) - u(cond3)));
            W_(cond4) = 0;

        case 'danska'
            q = c * dSigma;   
            cond = abs(d) <= q;
            W_(cond) = 1;
            W_(~cond) = exp(-(d(~cond) ./ q(~cond)).^2);
    end
    
 
    W_(isinf(W_)) = 1e9;
    W_(W_ > 1e9) = 1e9;
    W = diag(W_);

end

