function results = darm(dataset1, dataset2, wFunctionName, a, b, c, alpha, maxIterDiff)
% DARM Naredi deformacijsko analizo z izbrano robustno metodo
%
% Inputs:
%    dataset1      : vhodni podatki prve terminske izmere (iz .daf datoteke)
%    dataset2      : vhodni podatki druge terminske izmere (iz .daf datoteke)
%    wFunctionName : ime utežne funkcije (*L1*, L1-L2, Lp, Huber, HuberMod, Fair,
%                    Cauchy, Welsch, Tukey, German-McClure, Hampel, danska)
%    a, b, c       : parametri nekaterih utežnih funkcij
%    alpha         : stopnja značilnosti testa (*0.5*)
%    maxIterDiff   : največja razlika dy ozirmoa dx med dvema iteracijskima korakoma
%
% Outputs:
%     results      : struktura z rezultati


    arguments
        dataset1
        dataset2
        wFunctionName string {mustBeMember(wFunctionName, {'L1', 'L1-L2', 'Lp', 'Huber', 'HuberMod', ...
            'Fair', 'Cauchy', 'Welsch', 'Tukey', 'German-McClure', 'Hampel', 'danska'})} = 'L1'
        a {mustBeNumeric, mustBeReal} = -1
        b {mustBeNumeric, mustBeReal} = -1
        c {mustBeNumeric, mustBeReal} = -1
        alpha double = 0.05
        maxIterDiff {mustBeNumeric, mustBeReal} = 1e-8;
    end
    
    % skupno število nadštevilnosti in povprečna referenčna varianca a-posteriori (splošna aritmetična sredina)
    redundancyNumber = dataset1.redundancyNumber + dataset2.redundancyNumber;
    s0s0Apost = (dataset1.redundancyNumber * dataset1.s0Apost^2 + dataset2.redundancyNumber * dataset2.s0Apost^2) / redundancyNumber;
    
    % nastavimo začetni vektor premikov d = [dy1, dx1, dy2, dx2, dy3, dx3 ...]'
    dy = [dataset2.points.yAdj] - [dataset1.points.yAdj];
    dx = [dataset2.points.xAdj] - [dataset1.points.xAdj];
    nPoints = dataset1.nPoints;
    pointNames = [dataset1.points.name]';
    d = zeros(2*nPoints, 1);    
    d(1:2:end-1, 1) = dy;
    d(2:2:end  , 1) = dx;
    
    % nastavimo začetno matriko kofaktorjev Qdd
    Pdd = dataset1.matrixN * pinv(dataset1.matrixN + dataset2.matrixN) * dataset2.matrixN;
    Qdd = pinv(Pdd);
    %Qdd = dataset1.matrixQ + dataset2.matrixQ;
    Sdd = s0s0Apost * Qdd;
    
    % nastavimo datumsko matriko notranjih vezi H
    H = getMatrixH([dataset1.points.yApprox;dataset1.points.xApprox]');
    
    % ----------------------------------------------------------------------------------------------------
    % ZAČETNE VREDNOSTI ITERACIJSKEGA POSTOPKA
    % Če imamo prosto mrežo, ni potrebno nastavljati začetnih vrednosti, saj velja W = I --> d_ii = d_i.
    % W je nadsavljena za prosto mrežo.
    % Če mreža ni izravnana kot prosta:
    % W = D*inv(D'*D)*D'
    % S = eye(2*nPts) - H * ((H' * W * H) \ H') * W;
    % d_i = S * d;
    % dDiff = inf;
    % ----------------------------------------------------------------------------------------------------
    
    % ----------------------------------------------------------------------------------------------------
    % ---------------------------------- STABILNOST PO KOMPONENTAH (PC)-----------------------------------
    % ----------------------------------------------------------------------------------------------------
    d_i = d;
    dDiff = inf;
    nIterPC = 0;
    while max(abs(dDiff)) > maxIterDiff
        W = getMatrixW(d_i, sqrt(diag(Sdd)), wFunctionName, a, b, c);
        
        S = eye(2*nPoints) - H * ((H' * W * H) \ H') * W;
        d_ii = S * d_i;
        dDiff = d_ii - d_i;
        d_i = d_ii;
        
        nIterPC = nIterPC + 1;
        if nIterPC > 100
            error('Funkcija %s: Več kot 100 iteracij - iteracijski postopek prilagajanja uteži ne konvergira.', mfilename)
        end
    end
    
    % končni vektor premikov [dy1, dx1, dy2, dx2, dy3, dx3 ...]'
    dPC = d_i;
    
    % dolžine končnih vektorjev premikov s = [norm([dy1, dx1]) norm([dy2, dx2]) norm([dy3, dx3]) ...]'
    sPC = zeros(nPoints,1);
    for i = 1 : nPoints
        sPC(i,1) = norm([dPC(2*i-1), dPC(2*i)]);
    end
    
    % končna matrika kofaktorjev Qdd
    QddPC = (eye(2*nPoints) - H * ((H'*W*H)\H') * W) * Qdd * S';


    % ----------------------------------------------------------------------------------------------------
    % ------------------------------------- STABILNOST PO TOCKAH (PP)-------------------------------------
    % ----------------------------------------------------------------------------------------------------
    dy = d_i(1:2:end);
    dx = d_i(2:2:end);
    dyVar =   Sdd(          1 : 4*nPoints+2 : end)';
    dxVar =   Sdd(2*nPoints+2 : 4*nPoints+2 : end)';
    dydxCov = Sdd(          2 : 4*nPoints+2 : end)';
    
    % vektor celotnih premikov (ne po koordinatnih komponentah)
    s = sqrt(dy.^2 + dx.^2);
    
    % vektor diagonalnih elementov kovariančne matrike Sss
    sVar = ( dy ./ s ).^2 .* dyVar  +  2 * dy .* dx ./ s.^2 .* dydxCov  +  ( dx ./ s ).^2 .* dxVar;    
    
    d_i = d;
    s_i = s;
    dDiff = 1;
    nIterPP = 0;
    while max(abs(dDiff)) > maxIterDiff
        W_ = diag(getMatrixW(s_i, sqrt(sVar), wFunctionName, a, b, c));
        W = diag(repelem(W_, 2));  % končna matrika W - ker delamo analizo po točkah, jo je potrebno "podvojiti"

        S = eye(2*nPoints) - H * ((H' * W * H) \ H') * W;
        d_ii = S * d_i;
        dDiff = d_ii - d_i;
        d_i = d_ii;

        nIterPP = nIterPP + 1;
        if nIterPP > 100
            error('Funkcija %s: Več kot 100 iteracij - iteracijski postopek prilagajanja uteži ne konvergira.', mfilename)
        end
    end
    
    % končni vektor premikov [dy1, dx1, dy2, dx2, dy3, dx3 ...]'
    dPP = d_i;

    % dolžine končnih vektorjev premikov s = [norm([dy1, dx1]) norm([dy2, dx2]) norm([dy3, dx3]) ...]'
    sPP = zeros(nPoints, 1);
    for i = 1 : nPoints
        sPP(i) = sqrt(dPP(2*i-1)^2 + dPP(2*i)^2);
    end
    
    % končna matrika kofaktorjev Qdd
    QddPP = (eye(2*nPoints) - H * ((H'*W*H)\H') * W) * Qdd * S';

    
    % statistično testiranje
    iAlpha = 1 - (1 - alpha)^(1/nPoints);
    tCritValuePC = finv(1 - iAlpha, 1, redundancyNumber);
    tCritValuePP = finv(1 - iAlpha, 2, redundancyNumber);
    
    tStatisticPC = zeros(nPoints, 2);
    tStatisticPP = zeros(nPoints, 1);
    for i = 1 : nPoints
        jy = 2 * i - 1;
        jx = 2 * i;
        tStatisticPC(i, 1) = dPC(jy)' * QddPC(jy, jy)^-1 * dPC(jy) / (1 * s0s0Apost);
        tStatisticPC(i, 2) = dPC(jx)' * QddPC(jx, jx)^-1 * dPC(jx) / (1 * s0s0Apost);
        tStatisticPP(i) = dPP(jy:jx)' * QddPP(jy:jx, jy:jx)^-1 * dPP(jy:jx) / (2 * s0s0Apost);
    end
    

    % končni rezultati
    results.pointNames = pointNames;
    results.wFunctionName = wFunctionName;
    results.alpha = alpha;
    results.aParameter = a;
    results.bParameter = b;
    results.cParameter = c;

    results.perComponents.d = dPC;
    results.perComponents.tStatistic = tStatisticPC;
    results.perComponents.tCritValue = tCritValuePC;
    results.perComponents.stableComponents = tStatisticPC < tCritValuePC;
    results.perComponents.s = sPC;
    results.perComponents.stablePoints = results.perComponents.stableComponents(:,1) & results.perComponents.stableComponents(:,2);
    results.perComponents.nIter = nIterPC;
    
    results.perPoints.d = dPP;
    results.perPoints.s = sPP;
    results.perPoints.tStatistic = tStatisticPP;
    results.perPoints.tCritValue = tCritValuePP;
    results.perPoints.stablePoints = tStatisticPP < tCritValuePP;
    results.perPoints.nIter = nIterPP;
end

