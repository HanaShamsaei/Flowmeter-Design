function [P1, P2] = simulation(D_mm, d_mm, U0)
    import com.comsol.model.*
    import com.comsol.model.util.*
    
    % --- Fixed geometry ----------------------------------------------------
    
    D  = D_mm;         % [mm] outer diameter  (kept constant)
    d  = d_mm;       % [mm] inner diameter  (β = 0.5)
    U0in = U0;

    if D == 100
        alpha = 0.015;
    elseif D == 150
        alpha = 0.01;
    end

    e  = alpha*D;      % [mm] orifice thickness

    E  = 2*alpha*D;      % [mm] plate thickness 

    
    info = { ...
        'tbl1','Upstream flow rate  [m^3/s]'; ...
        'tbl2','Downstream flow rate [m^3/s]'; ...
        'tbl3','Avg axial velocity up [m/s]'; ...
        'tbl4','Avg pressure up [Pa]'; ...
        'tbl5','Avg axial velocity down [m/s]'; ...
        'tbl6','Avg pressure down [Pa]' };
    
    n = size(info,1);
    

        try
            % ----- build a fresh model --------------------------------------
            model = ModelUtil.create('model');
        
            model.modelPath('E:\C\University\Semesters\Semester 8\Bachelor''s Project\comsol\my simulations');
            
            model.component.create('comp1', true);
            
            model.component('comp1').geom.create('geom1', 2);
            model.component('comp1').geom('geom1').axisymmetric(true);
            
            model.component('comp1').mesh.create('mesh1');
            
            model.component('comp1').physics.create('spf', 'TurbulentFlowSST', 'geom1');
            
            model.study.create('std1');
            model.study('std1').create('wdi', 'WallDistanceInitialization');
            model.study('std1').feature('wdi').set('solnum', 'auto');
            model.study('std1').feature('wdi').set('notsolnum', 'auto');
            model.study('std1').feature('wdi').set('outputmap', {});
            model.study('std1').feature('wdi').set('ngenAUX', '1');
            model.study('std1').feature('wdi').set('goalngenAUX', '1');
            model.study('std1').feature('wdi').set('ngenAUX', '1');
            model.study('std1').feature('wdi').set('goalngenAUX', '1');
            model.study('std1').feature('wdi').setSolveFor('/physics/spf', true);
            model.study('std1').create('stat', 'Stationary');
            model.study('std1').feature('stat').set('solnum', 'auto');
            model.study('std1').feature('stat').set('notsolnum', 'auto');
            model.study('std1').feature('stat').set('outputmap', {});
            model.study('std1').feature('stat').set('ngenAUX', '1');
            model.study('std1').feature('stat').set('goalngenAUX', '1');
            model.study('std1').feature('stat').set('ngenAUX', '1');
            model.study('std1').feature('stat').set('goalngenAUX', '1');
            model.study('std1').feature('stat').setSolveFor('/physics/spf', true);
    
    
            model.component('comp1').geom('geom1').lengthUnit('mm');
            model.component('comp1').geom('geom1').create('r1', 'Rectangle');
            model.component('comp1').geom('geom1').feature('r1').set('size', [D/2 5000]);
            model.component('comp1').geom('geom1').run('r1');
            model.component('comp1').geom('geom1').create('r2', 'Rectangle');
            model.component('comp1').geom('geom1').feature('r2').set('size', [(D-d)/2 E]);
            model.component('comp1').geom('geom1').feature('r2').set('pos', [d/2 2500-E/2]);
            model.component('comp1').geom('geom1').runPre('fin');
            model.component('comp1').geom('geom1').create('cha1', 'Chamfer');
            model.component('comp1').geom('geom1').feature('cha1').selection('pointinsketch').set('r2', 4);
            model.component('comp1').geom('geom1').feature('cha1').set('dist', E-e);
            model.component('comp1').geom('geom1').runPre('fin');
            model.component('comp1').geom('geom1').create('dif1', 'Difference');
            model.component('comp1').geom('geom1').feature('dif1').selection('input').set({'r1'});
            model.component('comp1').geom('geom1').feature('dif1').selection('input2').set({'cha1'});
            model.component('comp1').geom('geom1').runPre('fin');
            model.component('comp1').geom('geom1').create('ls1', 'LineSegment');
            model.component('comp1').geom('geom1').feature('ls1').set('specify1', 'coord');
            model.component('comp1').geom('geom1').feature('ls1').set('coord1', [0 2500-D]);
            model.component('comp1').geom('geom1').feature('ls1').set('specify2', 'coord');
            model.component('comp1').geom('geom1').feature('ls1').set('coord2', [D/2 2500-D]);
            model.component('comp1').geom('geom1').runPre('fin');
            model.component('comp1').geom('geom1').create('ls2', 'LineSegment');
            model.component('comp1').geom('geom1').feature('ls2').set('specify1', 'coord');
            model.component('comp1').geom('geom1').feature('ls2').set('coord1', [0 2500+D/2]);
            model.component('comp1').geom('geom1').feature('ls2').set('specify2', 'coord');
            model.component('comp1').geom('geom1').feature('ls2').set('coord2', [D/2 2500+D/2]);
            model.component('comp1').geom('geom1').runPre('fin');
            model.component('comp1').geom('geom1').run;
    
    
            model.component('comp1').material.create('mat1', 'Common');
            model.component('comp1').material('mat1').propertyGroup('def').func.create('eta', 'Piecewise');
            model.component('comp1').material('mat1').propertyGroup('def').func.create('Cp', 'Piecewise');
            model.component('comp1').material('mat1').propertyGroup('def').func.create('rho', 'Piecewise');
            model.component('comp1').material('mat1').propertyGroup('def').func.create('k', 'Piecewise');
            model.component('comp1').material('mat1').propertyGroup('def').func.create('cs', 'Interpolation');
            model.component('comp1').material('mat1').propertyGroup('def').func.create('an1', 'Analytic');
            model.component('comp1').material('mat1').propertyGroup('def').func.create('an2', 'Analytic');
            model.component('comp1').material('mat1').propertyGroup('def').func.create('an3', 'Analytic');
            model.component('comp1').material('mat1').label('Water');
            model.component('comp1').material('mat1').set('family', 'water');
            model.component('comp1').material('mat1').propertyGroup('def').func('eta').set('arg', 'T');
            model.component('comp1').material('mat1').propertyGroup('def').func('eta').set('pieces', {'273.15' '413.15' '1.3799566804-0.021224019151*T^1+1.3604562827E-4*T^2-4.6454090319E-7*T^3+8.9042735735E-10*T^4-9.0790692686E-13*T^5+3.8457331488E-16*T^6'; '413.15' '553.75' '0.00401235783-2.10746715E-5*T^1+3.85772275E-8*T^2-2.39730284E-11*T^3'});
            model.component('comp1').material('mat1').propertyGroup('def').func('eta').set('argunit', 'K');
            model.component('comp1').material('mat1').propertyGroup('def').func('eta').set('fununit', 'Pa*s');
            model.component('comp1').material('mat1').propertyGroup('def').func('Cp').set('arg', 'T');
            model.component('comp1').material('mat1').propertyGroup('def').func('Cp').set('pieces', {'273.15' '553.75' '12010.1471-80.4072879*T^1+0.309866854*T^2-5.38186884E-4*T^3+3.62536437E-7*T^4'});
            model.component('comp1').material('mat1').propertyGroup('def').func('Cp').set('argunit', 'K');
            model.component('comp1').material('mat1').propertyGroup('def').func('Cp').set('fununit', 'J/(kg*K)');
            model.component('comp1').material('mat1').propertyGroup('def').func('rho').set('arg', 'T');
            model.component('comp1').material('mat1').propertyGroup('def').func('rho').set('smooth', 'contd1');
            model.component('comp1').material('mat1').propertyGroup('def').func('rho').set('pieces', {'273.15' '293.15' '0.000063092789034*T^3-0.060367639882855*T^2+18.9229382407066*T-950.704055329848'; '293.15' '373.15' '0.000010335053319*T^3-0.013395065634452*T^2+4.969288832655160*T+432.257114008512'});
            model.component('comp1').material('mat1').propertyGroup('def').func('rho').set('argunit', 'K');
            model.component('comp1').material('mat1').propertyGroup('def').func('rho').set('fununit', 'kg/m^3');
            model.component('comp1').material('mat1').propertyGroup('def').func('k').set('arg', 'T');
            model.component('comp1').material('mat1').propertyGroup('def').func('k').set('pieces', {'273.15' '553.75' '-0.869083936+0.00894880345*T^1-1.58366345E-5*T^2+7.97543259E-9*T^3'});
            model.component('comp1').material('mat1').propertyGroup('def').func('k').set('argunit', 'K');
            model.component('comp1').material('mat1').propertyGroup('def').func('k').set('fununit', 'W/(m*K)');
            model.component('comp1').material('mat1').propertyGroup('def').func('cs').set('table', {'273' '1403';  ...
            '278' '1427';  ...
            '283' '1447';  ...
            '293' '1481';  ...
            '303' '1507';  ...
            '313' '1526';  ...
            '323' '1541';  ...
            '333' '1552';  ...
            '343' '1555';  ...
            '353' '1555';  ...
            '363' '1550';  ...
            '373' '1543'});
            model.component('comp1').material('mat1').propertyGroup('def').func('cs').set('interp', 'piecewisecubic');
            model.component('comp1').material('mat1').propertyGroup('def').func('cs').set('fununit', {'m/s'});
            model.component('comp1').material('mat1').propertyGroup('def').func('cs').set('argunit', {'K'});
            model.component('comp1').material('mat1').propertyGroup('def').func('an1').label('Analytic ');
            model.component('comp1').material('mat1').propertyGroup('def').func('an1').set('funcname', 'alpha_p');
            model.component('comp1').material('mat1').propertyGroup('def').func('an1').set('expr', '-1/rho(T)*d(rho(T),T)');
            model.component('comp1').material('mat1').propertyGroup('def').func('an1').set('args', {'T'});
            model.component('comp1').material('mat1').propertyGroup('def').func('an1').set('fununit', '1/K');
            model.component('comp1').material('mat1').propertyGroup('def').func('an1').set('argunit', {'K'});
            model.component('comp1').material('mat1').propertyGroup('def').func('an1').set('plotfixedvalue', {'273.15'});
            model.component('comp1').material('mat1').propertyGroup('def').func('an1').set('plotargs', {'T' '273.15' '373.15'});
            model.component('comp1').material('mat1').propertyGroup('def').func('an2').set('funcname', 'gamma_w');
            model.component('comp1').material('mat1').propertyGroup('def').func('an2').set('expr', '1+(T/Cp(T))*(alpha_p(T)*cs(T))^2');
            model.component('comp1').material('mat1').propertyGroup('def').func('an2').set('args', {'T'});
            model.component('comp1').material('mat1').propertyGroup('def').func('an2').set('fununit', '1');
            model.component('comp1').material('mat1').propertyGroup('def').func('an2').set('argunit', {'K'});
            model.component('comp1').material('mat1').propertyGroup('def').func('an2').set('plotfixedvalue', {'273.15'});
            model.component('comp1').material('mat1').propertyGroup('def').func('an2').set('plotargs', {'T' '273.15' '373.15'});
            model.component('comp1').material('mat1').propertyGroup('def').func('an3').set('funcname', 'muB');
            model.component('comp1').material('mat1').propertyGroup('def').func('an3').set('expr', '2.79*eta(T)');
            model.component('comp1').material('mat1').propertyGroup('def').func('an3').set('args', {'T'});
            model.component('comp1').material('mat1').propertyGroup('def').func('an3').set('fununit', 'Pa*s');
            model.component('comp1').material('mat1').propertyGroup('def').func('an3').set('argunit', {'K'});
            model.component('comp1').material('mat1').propertyGroup('def').func('an3').set('plotfixedvalue', {'273.15'});
            model.component('comp1').material('mat1').propertyGroup('def').func('an3').set('plotargs', {'T' '273.15' '553.75'});
            model.component('comp1').material('mat1').propertyGroup('def').set('thermalexpansioncoefficient', '');
            model.component('comp1').material('mat1').propertyGroup('def').set('bulkviscosity', '');
            model.component('comp1').material('mat1').propertyGroup('def').set('thermalexpansioncoefficient', {'alpha_p(T)' '0' '0' '0' 'alpha_p(T)' '0' '0' '0' 'alpha_p(T)'});
            model.component('comp1').material('mat1').propertyGroup('def').set('bulkviscosity', 'muB(T)');
            model.component('comp1').material('mat1').propertyGroup('def').set('dynamicviscosity', 'eta(T)');
            model.component('comp1').material('mat1').propertyGroup('def').set('ratioofspecificheat', 'gamma_w(T)');
            model.component('comp1').material('mat1').propertyGroup('def').set('electricconductivity', {'5.5e-6[S/m]' '0' '0' '0' '5.5e-6[S/m]' '0' '0' '0' '5.5e-6[S/m]'});
            model.component('comp1').material('mat1').propertyGroup('def').set('heatcapacity', 'Cp(T)');
            model.component('comp1').material('mat1').propertyGroup('def').set('density', 'rho(T)');
            model.component('comp1').material('mat1').propertyGroup('def').set('thermalconductivity', {'k(T)' '0' '0' '0' 'k(T)' '0' '0' '0' 'k(T)'});
            model.component('comp1').material('mat1').propertyGroup('def').set('soundspeed', 'cs(T)');
            model.component('comp1').material('mat1').propertyGroup('def').addInput('temperature');
    
            model.component('comp1').physics('spf').create('inl1', 'InletBoundary', 1);
            model.component('comp1').physics('spf').feature('inl1').selection.set([2]);
            model.component('comp1').physics('spf').feature('inl1').set('U0in', U0in);
            model.component('comp1').physics('spf').create('out1', 'OutletBoundary', 1);
            model.component('comp1').physics('spf').feature('out1').selection.set([7]);
    
            model.component('comp1').mesh('mesh1').run;
    
            model.sol.create('sol1');
            
            model.component('comp1').mesh('mesh1').stat.selection.geom(2);
            model.component('comp1').mesh('mesh1').stat.selection.set([1 2 3]);
            model.component('comp1').mesh('mesh1').stat.selection.geom(2);
            model.component('comp1').mesh('mesh1').stat.selection.set([1 2 3]);
            model.component('comp1').mesh('mesh1').stat.selection.geom(2);
            model.component('comp1').mesh('mesh1').stat.selection.set([1 2 3]);
            
            model.sol('sol1').study('std1');
            model.sol('sol1').create('st1', 'StudyStep');
            model.sol('sol1').feature('st1').set('study', 'std1');
            model.sol('sol1').feature('st1').set('studystep', 'wdi');
            model.sol('sol1').create('v1', 'Variables');
            model.sol('sol1').feature('v1').set('control', 'wdi');
            model.sol('sol1').create('s1', 'Stationary');
            model.sol('sol1').feature('s1').set('stol', 1.0E-6);
            model.sol('sol1').feature('s1').feature('aDef').set('cachepattern', true);
            model.sol('sol1').feature('s1').create('fc1', 'FullyCoupled');
            model.sol('sol1').feature('s1').feature('fc1').set('dtech', 'auto');
            model.sol('sol1').feature('s1').feature('fc1').set('initstep', 0.01);
            model.sol('sol1').feature('s1').feature('fc1').set('minstep', 1.0E-6);
            model.sol('sol1').feature('s1').feature('fc1').set('maxiter', 50);
            model.sol('sol1').feature('s1').create('d1', 'Direct');
            model.sol('sol1').feature('s1').feature('d1').set('linsolver', 'pardiso');
            model.sol('sol1').feature('s1').feature('d1').set('pivotperturb', 1.0E-13);
            model.sol('sol1').feature('s1').feature('d1').label('Direct, wall distance (spf)');
            model.sol('sol1').feature('s1').create('i1', 'Iterative');
            model.sol('sol1').feature('s1').feature('i1').set('linsolver', 'gmres');
            model.sol('sol1').feature('s1').feature('i1').set('prefuntype', 'left');
            model.sol('sol1').feature('s1').feature('i1').set('itrestart', 50);
            model.sol('sol1').feature('s1').feature('i1').set('rhob', 400);
            model.sol('sol1').feature('s1').feature('i1').set('maxlinit', 1000);
            model.sol('sol1').feature('s1').feature('i1').set('nlinnormuse', 'on');
            model.sol('sol1').feature('s1').feature('i1').label('AMG, wall distance (spf)');
            model.sol('sol1').feature('s1').feature('i1').create('mg1', 'Multigrid');
            model.sol('sol1').feature('s1').feature('i1').feature('mg1').set('prefun', 'saamg');
            model.sol('sol1').feature('s1').feature('i1').feature('mg1').set('mgcycle', 'v');
            model.sol('sol1').feature('s1').feature('i1').feature('mg1').set('maxcoarsedof', 50000);
            model.sol('sol1').feature('s1').feature('i1').feature('mg1').set('strconn', 0.01);
            model.sol('sol1').feature('s1').feature('i1').feature('mg1').set('nullspace', 'constant');
            model.sol('sol1').feature('s1').feature('i1').feature('mg1').set('usesmooth', false);
            model.sol('sol1').feature('s1').feature('i1').feature('mg1').set('saamgcompwise', true);
            model.sol('sol1').feature('s1').feature('i1').feature('mg1').set('loweramg', true);
            model.sol('sol1').feature('s1').feature('i1').feature('mg1').set('compactaggregation', false);
            model.sol('sol1').feature('s1').feature('i1').feature('mg1').feature('pr').create('sl1', 'SORLine');
            model.sol('sol1').feature('s1').feature('i1').feature('mg1').feature('pr').feature('sl1').set('linesweeptype', 'ssor');
            model.sol('sol1').feature('s1').feature('i1').feature('mg1').feature('pr').feature('sl1').set('iter', 1);
            model.sol('sol1').feature('s1').feature('i1').feature('mg1').feature('pr').feature('sl1').set('linerelax', 0.7);
            model.sol('sol1').feature('s1').feature('i1').feature('mg1').feature('pr').feature('sl1').set('linealgorithm', 'mesh');
            model.sol('sol1').feature('s1').feature('i1').feature('mg1').feature('pr').feature('sl1').set('linemethod', 'uncoupled');
            model.sol('sol1').feature('s1').feature('i1').feature('mg1').feature('pr').feature('sl1').set('seconditer', 1);
            model.sol('sol1').feature('s1').feature('i1').feature('mg1').feature('pr').feature('sl1').set('relax', 0.5);
            model.sol('sol1').feature('s1').feature('i1').feature('mg1').feature('po').create('sl1', 'SORLine');
            model.sol('sol1').feature('s1').feature('i1').feature('mg1').feature('po').feature('sl1').set('linesweeptype', 'ssor');
            model.sol('sol1').feature('s1').feature('i1').feature('mg1').feature('po').feature('sl1').set('iter', 1);
            model.sol('sol1').feature('s1').feature('i1').feature('mg1').feature('po').feature('sl1').set('linerelax', 0.7);
            model.sol('sol1').feature('s1').feature('i1').feature('mg1').feature('po').feature('sl1').set('linealgorithm', 'mesh');
            model.sol('sol1').feature('s1').feature('i1').feature('mg1').feature('po').feature('sl1').set('linemethod', 'uncoupled');
            model.sol('sol1').feature('s1').feature('i1').feature('mg1').feature('po').feature('sl1').set('seconditer', 1);
            model.sol('sol1').feature('s1').feature('i1').feature('mg1').feature('po').feature('sl1').set('relax', 0.5);
            model.sol('sol1').feature('s1').feature('i1').feature('mg1').feature('cs').create('d1', 'Direct');
            model.sol('sol1').feature('s1').feature('i1').feature('mg1').feature('cs').feature('d1').set('linsolver', 'pardiso');
            model.sol('sol1').feature('s1').feature('i1').feature('mg1').feature('cs').feature('d1').set('pivotperturb', 1.0E-13);
            model.sol('sol1').feature('s1').feature('fc1').set('linsolver', 'd1');
            model.sol('sol1').feature('s1').feature('fc1').set('dtech', 'auto');
            model.sol('sol1').feature('s1').feature('fc1').set('initstep', 0.01);
            model.sol('sol1').feature('s1').feature('fc1').set('minstep', 1.0E-6);
            model.sol('sol1').feature('s1').feature('fc1').set('maxiter', 50);
            model.sol('sol1').feature('s1').feature.remove('fcDef');
            model.sol('sol1').create('su1', 'StoreSolution');
            model.sol('sol1').create('st2', 'StudyStep');
            model.sol('sol1').feature('st2').set('study', 'std1');
            model.sol('sol1').feature('st2').set('studystep', 'stat');
            model.sol('sol1').create('v2', 'Variables');
            model.sol('sol1').feature('v2').set('initmethod', 'sol');
            model.sol('sol1').feature('v2').set('initsol', 'sol1');
            model.sol('sol1').feature('v2').set('notsolmethod', 'sol');
            model.sol('sol1').feature('v2').set('notsol', 'sol1');
            model.sol('sol1').feature('v2').set('control', 'stat');
            model.sol('sol1').create('s2', 'Stationary');
            model.sol('sol1').feature('s2').feature('aDef').set('cachepattern', true);
            model.sol('sol1').feature('s2').create('se1', 'Segregated');
            model.sol('sol1').feature('s2').feature('se1').feature.remove('ssDef');
            model.sol('sol1').feature('s2').feature('se1').create('ss1', 'SegregatedStep');
            model.sol('sol1').feature('s2').feature('se1').feature('ss1').set('segvar', {'comp1_p' 'comp1_u'});
            model.sol('sol1').feature('s2').feature('se1').feature('ss1').set('subdamp', 0.5);
            model.sol('sol1').feature('s2').create('d1', 'Direct');
            model.sol('sol1').feature('s2').feature('d1').set('linsolver', 'pardiso');
            model.sol('sol1').feature('s2').feature('d1').set('pivotperturb', 1.0E-13);
            model.sol('sol1').feature('s2').feature('d1').label('Direct, fluid flow variables (spf)');
            model.sol('sol1').feature('s2').feature('se1').feature('ss1').set('linsolver', 'd1');
            model.sol('sol1').feature('s2').feature('se1').feature('ss1').label('Velocity u, Pressure p');
            model.sol('sol1').feature('s2').feature('se1').create('ss2', 'SegregatedStep');
            model.sol('sol1').feature('s2').feature('se1').feature('ss2').set('segvar', {'comp1_k' 'comp1_om'});
            model.sol('sol1').feature('s2').feature('se1').feature('ss2').set('subdamp', 0.45);
            model.sol('sol1').feature('s2').feature('se1').feature('ss2').set('subiter', 3);
            model.sol('sol1').feature('s2').feature('se1').feature('ss2').set('subtermconst', 'itertol');
            model.sol('sol1').feature('s2').feature('se1').feature('ss2').set('subntolfact', 1);
            model.sol('sol1').feature('s2').create('d2', 'Direct');
            model.sol('sol1').feature('s2').feature('d2').set('linsolver', 'pardiso');
            model.sol('sol1').feature('s2').feature('d2').set('pivotperturb', 1.0E-13);
            model.sol('sol1').feature('s2').feature('d2').label('Direct, turbulence variables (spf)');
            model.sol('sol1').feature('s2').feature('se1').feature('ss2').set('linsolver', 'd2');
            model.sol('sol1').feature('s2').feature('se1').feature('ss2').label('Turbulence Variables');
            model.sol('sol1').feature('s2').feature('se1').set('segstabacc', 'segcflcmp');
            model.sol('sol1').feature('s2').feature('se1').set('subinitcfl', 2);
            model.sol('sol1').feature('s2').feature('se1').set('submincfl', 10000);
            model.sol('sol1').feature('s2').feature('se1').set('subkppid', 0.65);
            model.sol('sol1').feature('s2').feature('se1').set('subkdpid', 0.05);
            model.sol('sol1').feature('s2').feature('se1').set('subkipid', 0.05);
            model.sol('sol1').feature('s2').feature('se1').set('subcfltol', 0.1);
            model.sol('sol1').feature('s2').feature('se1').set('segcflaa', true);
            model.sol('sol1').feature('s2').feature('se1').set('segcflaacfl', 9000);
            model.sol('sol1').feature('s2').feature('se1').set('segcflaafact', 1);
            model.sol('sol1').feature('s2').feature('se1').set('maxsegiter', 300);
            model.sol('sol1').feature('s2').feature('se1').create('ll1', 'LowerLimit');
            model.sol('sol1').feature('s2').feature('se1').feature('ll1').set('lowerlimit', 'comp1.k 0 comp1.om 0 ');
            model.sol('sol1').feature('s2').create('i1', 'Iterative');
            model.sol('sol1').feature('s2').feature('i1').set('linsolver', 'gmres');
            model.sol('sol1').feature('s2').feature('i1').set('prefuntype', 'left');
            model.sol('sol1').feature('s2').feature('i1').set('itrestart', 50);
            model.sol('sol1').feature('s2').feature('i1').set('rhob', 20);
            model.sol('sol1').feature('s2').feature('i1').set('maxlinit', 1000);
            model.sol('sol1').feature('s2').feature('i1').set('nlinnormuse', 'on');
            model.sol('sol1').feature('s2').feature('i1').label('AMG, fluid flow variables (spf)');
            model.sol('sol1').feature('s2').feature('i1').create('mg1', 'Multigrid');
            model.sol('sol1').feature('s2').feature('i1').feature('mg1').set('prefun', 'saamg');
            model.sol('sol1').feature('s2').feature('i1').feature('mg1').set('mgcycle', 'v');
            model.sol('sol1').feature('s2').feature('i1').feature('mg1').set('maxcoarsedof', 80000);
            model.sol('sol1').feature('s2').feature('i1').feature('mg1').set('strconn', 0.02);
            model.sol('sol1').feature('s2').feature('i1').feature('mg1').set('nullspace', 'constant');
            model.sol('sol1').feature('s2').feature('i1').feature('mg1').set('usesmooth', false);
            model.sol('sol1').feature('s2').feature('i1').feature('mg1').set('saamgcompwise', true);
            model.sol('sol1').feature('s2').feature('i1').feature('mg1').set('loweramg', true);
            model.sol('sol1').feature('s2').feature('i1').feature('mg1').set('compactaggregation', false);
            model.sol('sol1').feature('s2').feature('i1').feature('mg1').feature('pr').create('sc1', 'SCGS');
            model.sol('sol1').feature('s2').feature('i1').feature('mg1').feature('pr').feature('sc1').set('linesweeptype', 'ssor');
            model.sol('sol1').feature('s2').feature('i1').feature('mg1').feature('pr').feature('sc1').set('iter', 0);
            model.sol('sol1').feature('s2').feature('i1').feature('mg1').feature('pr').feature('sc1').set('scgsrelax', 0.7);
            model.sol('sol1').feature('s2').feature('i1').feature('mg1').feature('pr').feature('sc1').set('scgsmethod', 'lines');
            model.sol('sol1').feature('s2').feature('i1').feature('mg1').feature('pr').feature('sc1').set('scgsvertexrelax', 0.7);
            model.sol('sol1').feature('s2').feature('i1').feature('mg1').feature('pr').feature('sc1').set('relax', 0.5);
            model.sol('sol1').feature('s2').feature('i1').feature('mg1').feature('pr').feature('sc1').set('scgssolv', 'stored');
            model.sol('sol1').feature('s2').feature('i1').feature('mg1').feature('pr').feature('sc1').set('approxscgs', true);
            model.sol('sol1').feature('s2').feature('i1').feature('mg1').feature('pr').feature('sc1').set('scgsdirectmaxsize', 1000);
            model.sol('sol1').feature('s2').feature('i1').feature('mg1').feature('po').create('sc1', 'SCGS');
            model.sol('sol1').feature('s2').feature('i1').feature('mg1').feature('po').feature('sc1').set('linesweeptype', 'ssor');
            model.sol('sol1').feature('s2').feature('i1').feature('mg1').feature('po').feature('sc1').set('iter', 1);
            model.sol('sol1').feature('s2').feature('i1').feature('mg1').feature('po').feature('sc1').set('scgsrelax', 0.7);
            model.sol('sol1').feature('s2').feature('i1').feature('mg1').feature('po').feature('sc1').set('scgsmethod', 'lines');
            model.sol('sol1').feature('s2').feature('i1').feature('mg1').feature('po').feature('sc1').set('scgsvertexrelax', 0.7);
            model.sol('sol1').feature('s2').feature('i1').feature('mg1').feature('po').feature('sc1').set('relax', 0.5);
            model.sol('sol1').feature('s2').feature('i1').feature('mg1').feature('po').feature('sc1').set('scgssolv', 'stored');
            model.sol('sol1').feature('s2').feature('i1').feature('mg1').feature('po').feature('sc1').set('approxscgs', true);
            model.sol('sol1').feature('s2').feature('i1').feature('mg1').feature('po').feature('sc1').set('scgsdirectmaxsize', 1000);
            model.sol('sol1').feature('s2').feature('i1').feature('mg1').feature('cs').create('d1', 'Direct');
            model.sol('sol1').feature('s2').feature('i1').feature('mg1').feature('cs').feature('d1').set('linsolver', 'pardiso');
            model.sol('sol1').feature('s2').feature('i1').feature('mg1').feature('cs').feature('d1').set('pivotperturb', 1.0E-13);
            model.sol('sol1').feature('s2').create('i2', 'Iterative');
            model.sol('sol1').feature('s2').feature('i2').set('linsolver', 'gmres');
            model.sol('sol1').feature('s2').feature('i2').set('prefuntype', 'left');
            model.sol('sol1').feature('s2').feature('i2').set('itrestart', 50);
            model.sol('sol1').feature('s2').feature('i2').set('rhob', 20);
            model.sol('sol1').feature('s2').feature('i2').set('maxlinit', 1000);
            model.sol('sol1').feature('s2').feature('i2').set('nlinnormuse', 'on');
            model.sol('sol1').feature('s2').feature('i2').label('AMG, turbulence variables (spf)');
            model.sol('sol1').feature('s2').feature('i2').create('mg1', 'Multigrid');
            model.sol('sol1').feature('s2').feature('i2').feature('mg1').set('prefun', 'saamg');
            model.sol('sol1').feature('s2').feature('i2').feature('mg1').set('mgcycle', 'v');
            model.sol('sol1').feature('s2').feature('i2').feature('mg1').set('maxcoarsedof', 50000);
            model.sol('sol1').feature('s2').feature('i2').feature('mg1').set('strconn', 0.01);
            model.sol('sol1').feature('s2').feature('i2').feature('mg1').set('nullspace', 'constant');
            model.sol('sol1').feature('s2').feature('i2').feature('mg1').set('usesmooth', false);
            model.sol('sol1').feature('s2').feature('i2').feature('mg1').set('saamgcompwise', true);
            model.sol('sol1').feature('s2').feature('i2').feature('mg1').set('loweramg', true);
            model.sol('sol1').feature('s2').feature('i2').feature('mg1').set('compactaggregation', false);
            model.sol('sol1').feature('s2').feature('i2').feature('mg1').feature('pr').create('sl1', 'SORLine');
            model.sol('sol1').feature('s2').feature('i2').feature('mg1').feature('pr').feature('sl1').set('linesweeptype', 'ssor');
            model.sol('sol1').feature('s2').feature('i2').feature('mg1').feature('pr').feature('sl1').set('iter', 0);
            model.sol('sol1').feature('s2').feature('i2').feature('mg1').feature('pr').feature('sl1').set('linerelax', 0.7);
            model.sol('sol1').feature('s2').feature('i2').feature('mg1').feature('pr').feature('sl1').set('linealgorithm', 'mesh');
            model.sol('sol1').feature('s2').feature('i2').feature('mg1').feature('pr').feature('sl1').set('linemethod', 'uncoupled');
            model.sol('sol1').feature('s2').feature('i2').feature('mg1').feature('pr').feature('sl1').set('seconditer', 1);
            model.sol('sol1').feature('s2').feature('i2').feature('mg1').feature('pr').feature('sl1').set('relax', 0.5);
            model.sol('sol1').feature('s2').feature('i2').feature('mg1').feature('po').create('sl1', 'SORLine');
            model.sol('sol1').feature('s2').feature('i2').feature('mg1').feature('po').feature('sl1').set('linesweeptype', 'ssor');
            model.sol('sol1').feature('s2').feature('i2').feature('mg1').feature('po').feature('sl1').set('iter', 1);
            model.sol('sol1').feature('s2').feature('i2').feature('mg1').feature('po').feature('sl1').set('linerelax', 0.7);
            model.sol('sol1').feature('s2').feature('i2').feature('mg1').feature('po').feature('sl1').set('linealgorithm', 'mesh');
            model.sol('sol1').feature('s2').feature('i2').feature('mg1').feature('po').feature('sl1').set('linemethod', 'uncoupled');
            model.sol('sol1').feature('s2').feature('i2').feature('mg1').feature('po').feature('sl1').set('seconditer', 1);
            model.sol('sol1').feature('s2').feature('i2').feature('mg1').feature('po').feature('sl1').set('relax', 0.5);
            model.sol('sol1').feature('s2').feature('i2').feature('mg1').feature('cs').create('d1', 'Direct');
            model.sol('sol1').feature('s2').feature('i2').feature('mg1').feature('cs').feature('d1').set('linsolver', 'pardiso');
            model.sol('sol1').feature('s2').feature('i2').feature('mg1').feature('cs').feature('d1').set('pivotperturb', 1.0E-13);
            model.sol('sol1').feature('s2').feature.remove('fcDef');
            model.sol('sol1').feature('v2').set('notsolnum', 'auto');
            model.sol('sol1').feature('v2').set('notsolvertype', 'solnum');
            model.sol('sol1').feature('v2').set('solnum', 'auto');
            model.sol('sol1').feature('v2').set('solvertype', 'solnum');
            model.sol('sol1').attach('std1');
            model.sol('sol1').runAll;
    
            model.result.dataset('dset1').set('geom', 'geom1');
            model.result.create('pg1', 'PlotGroup2D');
            model.result('pg1').label('Velocity (spf)');
            model.result('pg1').set('dataisaxisym', 'off');
            model.result('pg1').set('frametype', 'spatial');
            model.result('pg1').set('defaultPlotID', 'ResultDefaults_SinglePhaseFlow/icom1/pdef1/pcond1/pg1');
            model.result('pg1').feature.create('surf1', 'Surface');
            model.result('pg1').feature('surf1').label('Surface');
            model.result('pg1').feature('surf1').set('showsolutionparams', 'on');
            model.result('pg1').feature('surf1').set('smooth', 'internal');
            model.result('pg1').feature('surf1').set('showsolutionparams', 'on');
            model.result('pg1').feature('surf1').set('data', 'parent');
            model.result.create('pg2', 'PlotGroup2D');
            model.result('pg2').label('Pressure (spf)');
            model.result('pg2').set('dataisaxisym', 'off');
            model.result('pg2').set('frametype', 'spatial');
            model.result('pg2').set('defaultPlotID', 'ResultDefaults_SinglePhaseFlow/icom1/pdef1/pcond1/pg2');
            model.result('pg2').feature.create('con1', 'Contour');
            model.result('pg2').feature('con1').label('Contour');
            model.result('pg2').feature('con1').set('showsolutionparams', 'on');
            model.result('pg2').feature('con1').set('expr', 'p');
            model.result('pg2').feature('con1').set('number', 40);
            model.result('pg2').feature('con1').set('levelrounding', false);
            model.result('pg2').feature('con1').set('smooth', 'internal');
            model.result('pg2').feature('con1').set('showsolutionparams', 'on');
            model.result('pg2').feature('con1').set('data', 'parent');
            model.result.dataset.create('rev1', 'Revolve2D');
            model.result.dataset('rev1').label('Revolution 2D');
            model.result.dataset('rev1').set('data', 'none');
            model.result.dataset('rev1').set('startangle', -90);
            model.result.dataset('rev1').set('revangle', 225);
            model.result.dataset('rev1').set('data', 'dset1');
            model.result.create('pg3', 'PlotGroup3D');
            model.result('pg3').label('Velocity, 3D (spf)');
            model.result('pg3').set('frametype', 'spatial');
            model.result('pg3').set('defaultPlotID', 'ResultDefaults_SinglePhaseFlow/icom1/pdef1/pcond1/pcond1/pg1');
            model.result('pg3').feature.create('surf1', 'Surface');
            model.result('pg3').feature('surf1').label('Surface');
            model.result('pg3').feature('surf1').set('showsolutionparams', 'on');
            model.result('pg3').feature('surf1').set('smooth', 'internal');
            model.result('pg3').feature('surf1').set('showsolutionparams', 'on');
            model.result('pg3').feature('surf1').set('data', 'parent');
            model.result.dataset.create('edg1', 'Edge2D');
            model.result.dataset('edg1').label('Exterior Walls');
            model.result.dataset('edg1').set('data', 'dset1');
            model.result.dataset('edg1').selection.geom('geom1', 1);
            model.result.dataset('edg1').selection.set([8 9 10 11 12 13 14 15]);
            model.result.create('pg4', 'PlotGroup2D');
            model.result('pg4').label('Wall Resolution (spf)');
            model.result('pg4').set('dataisaxisym', 'off');
            model.result('pg4').set('frametype', 'spatial');
            model.result('pg4').set('defaultPlotID', 'ResultDefaults_SinglePhaseFlow/icom1/pdef1/pcond1/pcond2/pcond4/pg2');
            model.result('pg4').feature.create('line1', 'Line');
            model.result('pg4').feature('line1').label('Wall Resolution');
            model.result('pg4').feature('line1').set('showsolutionparams', 'on');
            model.result('pg4').feature('line1').set('expr', 'spf.Delta_wPlus');
            model.result('pg4').feature('line1').set('linetype', 'tube');
            model.result('pg4').feature('line1').set('smooth', 'internal');
            model.result('pg4').feature('line1').set('showsolutionparams', 'on');
            model.result('pg4').feature('line1').set('data', 'parent');
            model.result('pg4').feature('line1').feature.create('hght1', 'Height');
            model.result('pg4').feature('line1').feature('hght1').label('Height Expression');
            model.result('pg4').feature('line1').feature('hght1').set('heightdata', 'expr');
            model.result('pg4').feature('line1').feature('hght1').set('expr', 'spf.WRHeightExpr');
            model.result('pg1').run;
            model.result.dataset.create('edg2', 'Edge2D');
            model.result.dataset('edg2').selection.set([1 3 5]);
            model.result.numerical.create('int1', 'IntLine');
            model.result.numerical('int1').set('intsurface', true);
            model.result.numerical('int1').selection.set([4]);
            model.result.numerical('int1').set('expr', {'w'});
            model.result.numerical('int1').set('descr', {'Velocity field, z-component'});
            model.result.numerical('int1').set('unit', {'m^3/s'});
            model.result.table.create('tbl1', 'Table');
            model.result.table('tbl1').comments('Line Integration 1');
            model.result.numerical('int1').set('table', 'tbl1');
            model.result.numerical('int1').setResult;
            model.result.numerical('int1').label('flow rate-upstream');
            model.result.numerical.create('int2', 'IntLine');
            model.result.numerical('int2').set('intsurface', true);
            model.result.numerical('int2').selection.set([6]);
            model.result.numerical('int2').set('expr', {'w'});
            model.result.numerical('int2').set('descr', {'Velocity field, z-component'});
            model.result.numerical('int2').set('unit', {'m^3/s'});
            model.result.table.create('tbl2', 'Table');
            model.result.table('tbl2').comments('Line Integration 2');
            model.result.numerical('int2').set('table', 'tbl2');
            model.result.numerical('int2').setResult;
            model.result.numerical('int2').label('flow rate-downstream');
            model.result.numerical.create('av1', 'AvLine');
            model.result.numerical('av1').set('intsurface', true);
            model.result.numerical('av1').selection.set([4]);
            model.result.numerical('av1').set('expr', {'w'});
            model.result.numerical('av1').set('descr', {'Velocity field, z-component'});
            model.result.numerical('av1').set('unit', {'m/s'});
            model.result.table.create('tbl3', 'Table');
            model.result.table('tbl3').comments('Line Average 1');
            model.result.numerical('av1').set('table', 'tbl3');
            model.result.numerical('av1').setResult;
            model.result.numerical.create('av2', 'AvLine');
            model.result.numerical('av2').set('intsurface', true);
            model.result.numerical('av2').selection.set([4]);
            model.result.numerical('av2').set('expr', {'p'});
            model.result.numerical('av2').set('descr', {'Pressure'});
            model.result.numerical('av2').set('unit', {'Pa'});
            model.result.table.create('tbl4', 'Table');
            model.result.table('tbl4').comments('Line Average 2');
            model.result.numerical('av2').set('table', 'tbl4');
            model.result.numerical('av2').setResult;
            model.result.numerical('av1').label('velocity-upstream');
            model.result.numerical('av2').label('pressure-upstream');
            model.result.numerical.create('av3', 'AvLine');
            model.result.numerical('av3').set('intsurface', true);
            model.result.numerical('av3').selection.set([6]);
            model.result.numerical('av3').set('expr', {'w'});
            model.result.numerical('av3').set('descr', {'Velocity field, z-component'});
            model.result.numerical('av3').set('unit', {'m/s'});
            model.result.table.create('tbl5', 'Table');
            model.result.table('tbl5').comments('Line Average 3');
            model.result.numerical('av3').set('table', 'tbl5');
            model.result.numerical('av3').setResult;
            model.result.numerical('av3').label('velocity-downstream');
            model.result.numerical.create('av4', 'AvLine');
            model.result.numerical('av4').set('intsurface', true);
            model.result.numerical('av4').label('pressure-downstream');
            model.result.numerical('av4').selection.set([6]);
            model.result.numerical('av4').set('expr', {'p'});
            model.result.numerical('av4').set('descr', {'Pressure'});
            model.result.numerical('av4').set('unit', {'Pa'});
            model.result.table.create('tbl6', 'Table');
            model.result.table('tbl6').comments('pressure-downstream');
            model.result.numerical('av4').set('table', 'tbl6');
            model.result.numerical('av4').setResult;
            model.result.create('pg5', 'PlotGroup1D');
            model.result('pg5').run;
            model.result('pg5').create('lngr1', 'LineGraph');
            model.result('pg5').feature('lngr1').set('markerpos', 'datapoints');
            model.result('pg5').feature('lngr1').set('linewidth', 'preference');
            model.result('pg5').feature('lngr1').set('data', 'edg2');
            model.result('pg5').feature('lngr1').set('expr', 'p');
            model.result('pg5').run;
            model.result('pg5').run;
            model.result('pg5').create('lngr2', 'LineGraph');
            model.result('pg5').feature('lngr2').set('markerpos', 'datapoints');
            model.result('pg5').feature('lngr2').set('linewidth', 'preference');
            model.result('pg5').feature('lngr2').set('data', 'edg2');
            model.result('pg5').run;
            model.result('pg5').run;
            model.result('pg5').run;
            model.result('pg5').run;
            
            out = model;
    
            values = nan(1,n);
            for k = 1:n
                TT = mphtable(model,info{k,1});
                values(k) = TT.data(end);
            end

            P1 = values(4);
            P2 = values(6);
            
            ModelUtil.remove('model');
    
        catch ME
            warning('Run failed for U0in = %.2f m/s: %s', U0in, ME.message);
            P1 = NaN; P2 = NaN;   % <- ensure outputs exist
            ModelUtil.remove('model');
        end
end