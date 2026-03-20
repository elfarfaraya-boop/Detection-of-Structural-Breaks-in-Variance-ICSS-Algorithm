proc import datafile="C:\Users\ayael\Downloads\USD_RUB Historical Data (1).csv"
    out=rub_raw dbms=csv replace;
    getnames=yes;
run;

data rub_clean;
    set rub_raw (keep=Date Price);
    price_clean = input(compress(Price, ','), best12.);
    date_fmt = input(Date, ddmmyy10.);
    format date_fmt date9.;
    keep date_fmt price_clean;
    rename date_fmt   = date
           price_clean = close;
run;

proc sort data=rub_clean; by date; run;

/*Calcul des rendements logarithmique*/

data rendements;
    set rub_clean;
    close_lag = lag(close);
    if _n_ > 1 and close_lag > 0 then
        r = log(close / close_lag);
    if _n_ > 1;
    keep date close r;
run;

/*Statistiques descriptives*/

proc means data=rendements mean std skewness kurtosis;
    var r;
    title "Statistiques descriptives — USD/RUB";
run;

/*Test de stationnarité*/

proc arima data=rendements;
    identify var=r stationarity=(adf=5);
    title "Test ADF — Rendements USD/RUB";
run; quit;

proc arima data=rendements;
    identify var=r nlag=24;
    estimate p=1 q=1 method=ml;
    output out=residus residual=resid;
run; quit;


/*Procédure ICSS*/
proc iml;
    load module=icss_test;

    use residus;
    read all var {resid} into a;
    close residus;

    n_total    = nrow(a);
    a2         = a ## 2;
    D_critique = 1.224;

    cp_fwd = j(50, 1, 0);
    n_fwd  = 0;
    L = 1;
    do while(L < n_total);
        tau = icss_test(a2, L, n_total);
        if tau = 0 then leave;
        stable = 0;
        do while(stable = 0);
            tau2 = icss_test(a2, L, tau);
            if tau2 = 0 | tau2 = tau then stable = 1;
            else tau = tau2;
        end;
        n_fwd = n_fwd + 1;
        cp_fwd[n_fwd] = tau;
        L = tau + 1;
    end;

    cp_bwd = j(50, 1, 0);
    n_bwd  = 0;
    R = n_total;
    do while(R > 1);
        tau = icss_test(a2, 1, R);
        if tau = 0 then leave;
        stable = 0;
        do while(stable = 0);
            tau2 = icss_test(a2, tau + 1, R);
            if tau2 = 0 | tau2 = tau then stable = 1;
            else tau = tau2;
        end;
        n_bwd = n_bwd + 1;
        cp_bwd[n_bwd] = tau;
        R = tau;
    end;

    all_cp = j(100, 1, 0);
    n_all  = 0;
    do i = 1 to n_fwd;
        n_all = n_all + 1;
        all_cp[n_all] = cp_fwd[i];
    end;
    do i = 1 to n_bwd;
        n_all = n_all + 1;
        all_cp[n_all] = cp_bwd[i];
    end;
    if n_all > 0 then do;
        all_cp = all_cp[1:n_all];
        call sort(all_cp, 1);
        unique = j(n_all, 1, 0);
        unique[1] = all_cp[1];
        n_unique = 1;
        do i = 2 to n_all;
            if all_cp[i] ^= all_cp[i - 1] then do;
                n_unique = n_unique + 1;
                unique[n_unique] = all_cp[i];
            end;
        end;
        all_cp = unique[1:n_unique];
        n_all  = n_unique;
    end;

    converged = 0;
    do iter = 1 to 20 while(converged = 0);
        new_cp   = j(50, 1, 0);
        n_new    = 0;
        n_bounds = n_all + 2;
        bounds   = j(n_bounds, 1, 0);
        bounds[1] = 0;
        do i = 1 to n_all;
            bounds[i + 1] = all_cp[i];
        end;
        bounds[n_bounds] = n_total;
        do i = 1 to n_bounds - 1;
            seg_L = bounds[i] + 1;
            seg_R = bounds[i + 1];
            tau = icss_test(a2, seg_L, seg_R);
            if tau > 0 then do;
                n_new = n_new + 1;
                new_cp[n_new] = tau;
            end;
        end;
        if n_new = 0 & n_all = 0 then converged = 1;
        else if n_new = n_all then do;
            new_sorted = new_cp[1:n_new];
            call sort(new_sorted, 1);
            if all(new_sorted = all_cp) then converged = 1;
            else do; all_cp = new_sorted; n_all = n_new; end;
        end;
        else if n_new > 0 then do;
            all_cp = new_cp[1:n_new];
            call sort(all_cp, 1);
            n_all = n_new;
        end;
        else n_all = 0;
    end;

    if n_all > 0 then do;
        result = all_cp[1:n_all];
        print result[label="Points de rupture détectés"];
        print n_all[label="Nombre de ruptures"];
    end;
    else print "Aucune rupture détectée";

    C   = cusum(a2);
    C_N = C[n_total];
    k   = T(1:n_total);
    D   = (C / C_N) - (k / n_total);

    graphmat = k || a || C || D;
    create graphdata from graphmat[colname={"t" "a_t" "C_k" "D_k"}];
    append from graphmat;
    close graphdata;

    if n_all > 0 then do;
        create cpdata from all_cp[colname={"cp"}];
        append from all_cp;
        close cpdata;
    end;
    else do;
        cp_vide = {.};
        create cpdata from cp_vide[colname={"cp"}];
        append from cp_vide;
        close cpdata;
    end;

quit;


proc sql noprint;
    select cp into :cplist separated by ' '
    from cpdata where cp ^= .;
    select count(cp) into :ncp trimmed
    from cpdata where cp ^= .;
quit;

/* Graphiques */
proc sgplot data=graphdata;
    series x=t y=a_t / lineattrs=(color=black thickness=1);
    %if &ncp > 0 %then %do;
        refline &cplist / axis=x lineattrs=(color=red pattern=dash thickness=2)
                          label=("Ruptures");
    %end;
    refline 0 / axis=y lineattrs=(color=gray pattern=dot);
    xaxis label="Observation";
    yaxis label="Résidus ARMA";
    title "Résidus ARMA — USD/RUB";
run;

proc sgplot data=graphdata;
    series x=t y=C_k / lineattrs=(color=blue thickness=2);
    %if &ncp > 0 %then %do;
        refline &cplist / axis=x lineattrs=(color=red pattern=dash thickness=2);
    %end;
    xaxis label="Observation";
    yaxis label="Ck";
    title "Somme cumulée des carrés (Ck) — USD/RUB";
run;

proc sgplot data=graphdata;
    series x=t y=D_k / lineattrs=(color=black thickness=2);
    refline 0 / axis=y lineattrs=(color=gray pattern=dash);
    %if &ncp > 0 %then %do;
        refline &cplist / axis=x lineattrs=(color=red pattern=dash thickness=2);
    %end;
    xaxis label="Observation";
    yaxis label="Dk";
    title "Statistique Dk — USD/RUB";
run;

/* Dates de ruptures */
data dates_ruptures;
    set cpdata;
    if cp ^= .;
    obs = cp;
    set rendements point=obs;
    format date date9.;
    keep cp date close;
run;

proc print data=dates_ruptures noobs;
    title "Dates des ruptures de variance — USD/RUB";
run;

