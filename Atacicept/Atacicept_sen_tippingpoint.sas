/*
本ページではSAP1.0版の主要な解析に対する感度解析である、Tipping point 解析を行っていきます。
流れは以下になります。

Step1：データの整形
主解析用の補完済みデータセットから、
・透析、腎移植、または死亡
・試験からの脱落
が生じた被験者の36週の値を欠測とする。

Step2：欠測の補完
Stap1で作成された各補完されたデータセットに対し、シフトパラメータのペア(k1,k2)を用いて、MNARの仮定のもとで補完を行う
* k1= -2.7, -2.5, …, 2.5, 2.7
* k2= -2.7, -2.5, …, 2.5, 2.7

Step3：統計解析
MMRMを実施し、Rubin's ruleで結果を統合

*/


/*
Step1前に、主要な解析と同様の補完済みデータを作成していきます。
*/

*データの読み込み;
proc import datafile="C:\Users\XXX\Desktop\Journal\Atacicept\analysis_data.csv" out=imported_data dbms=csv replace;/*適宜パスの設定*/
     guessingrows=max;
run;

*NAが原因でOutcome_obsが文字列になっているため数値型に変換;
data data_n;
    set imported_data;
    temp_char = Outcome_obs;
    if temp_char = "NA" then temp_char = " ";
    Outcome_num = input(temp_char, ?? best12.);
    drop Outcome_obs temp_char;
    rename Outcome_num = Outcome_obs;
run;

*IDとグループでソート;
proc sort data=data_n out = data_sort;
    by PatientID trt;
run;

*PROC MIのためデータを横持ちに変換;
proc transpose data=data_sort out=wide_data prefix=Week;
    by PatientID trt missing_reason;
    id Week;
    var Outcome_obs;
run;

proc sort data = wide_data out = wide_data_sort ;
	by trt;
run;

/*計算の都合上補完回数やmcmcの設定を減らしています*/
proc mi data=wide_data_sort seed=123 nimpute=30 out=mono;
  by trt;
  var Week0 WEEK12 WEEK24 WEEK36;
  mcmc impute=monotone nbiter=500 niter=500; 
run;

*補完後データの確認;
proc mi data=mono nimpute=0;
	by trt;
    var Week0 Week12 Week24 Week36;
run;


*出力データをみると、両群ともに間欠的な欠測が補完できています。;
*本来はburn in期間や反復回数などが十分か確認すべきところですが、今回は省略します。;

proc sort data = mono out = mono_sort;
	by _imputation_ ;
run;

ods select none;
proc mi data=mono_sort seed=123 nimpute=1 out=outmi;
	class trt;
	by _imputation_ ;
	monotone reg(/details);
	mnar model (WEEK12 / modelobs= (trt= "Placebo"));*プラセボ群をreferenceに指定;
	mnar model (WEEK24 / modelobs= (trt= "Placebo"));
	mnar model (WEEK36 / modelobs= (trt= "Placebo"));
	var  WEEK0 WEEK12 WEEK24 WEEK36;
run;

/*
データの構造が分かりづらいのですが、missing reasonが「Dialysis/Transplant/Death」あるいは「Dropout」となっている症例が、
・透析、腎移植、または死亡
・試験からの脱落
となります(outmiは補完済みのため、補完値が入っています)。
Wide_data_sortなども適宜ご確認ください。
＊間欠的な欠測が生じている症例はmissug reasonが「Observed」となっています。
ここまでが主要な解析と同様の補完になります。
*/

/*
Step1：データの整形
主解析用の補完済みデータセットから、
・透析、腎移植、または死亡
・試験からの脱落
が生じた被験者の36週の値を欠測としていきます。
*/

data mi_tip_temp;
	set outmi;
	if missing_reason = "Dialysis/Transplant/Death" or missing_reason = "Dropout" then Week36_tip = .; else Week36_tip = Week36;
run;


/*Step2：欠測の補完
Stap1で作成されたデータセットに対し、シフトパラメータのペア(k1,k2)を用いて、MNARの仮定のもとで補完をしていきます(デルタ調整法)。
論文では-2.7から2.7まで0.2刻みなのですが、計算時間の都合上-3から3までの0.5刻みで設定。
* k1= -3.0, ... , 3.0
* k2= -3.0, ... , 3.0
*/

options nosource nonotes;
%let klist = -3.0 -2.5 -2.0 -1.5 -1.0 -0.5 0 0.5 1.0 1.5 2.0 2.5 3.0;

%macro tipping_mi_analysis;

%local i j k1 k2;

/* 以前の結果が残らないよう、最終出力用データセットをクリア */
proc datasets lib=work nolist;
    delete all_estimates;
quit;

%do i = 1 %to %sysfunc(countw(&klist, %str( )));
    %let k1 = %scan(&klist, &i, %str( ));/*区切り文字を空白のみに指定(これをしないと、3.0が3と0のように分割して認識されました)*/

    %do j = 1 %to %sysfunc(countw(&klist, %str( )));
       %let k2 = %scan(&klist, &j, %str( ));

        /* 1. 多重代入 (PROC MI) */
        proc mi data=mi_tip_temp nimpute=1 seed=%eval(1234 + &i*100 + &j) out=mi_tip_temp_out noprint;
            by _Imputation_;
            class trt;
            var WEEK0 WEEK12 WEEK24 WEEK36_tip;
            monotone reg;
            mnar adjust (WEEK36_tip / shift=&k1 adjustobs=(trt="Placebo"));
            mnar adjust (WEEK36_tip / shift=&k2 adjustobs=(trt="Atacicept"));
        run;

        /* 2. MMRMのため縦型に変換  */
        proc sort data=mi_tip_temp_out;
            by _imputation_ PatientID trt;
        run;

        proc transpose data=mi_tip_temp_out out=outmi_long name=VISIT;
            by _imputation_ PatientID trt;
            var WEEK0 WEEK12 WEEK24 WEEK36_tip; 
        run;

        /* 3. 解析  */
        proc mixed data=outmi_long(where=(VISIT ne "Week0")) method=reml;
            by _imputation_;
            class trt(ref = "Placebo") PatientID VISIT;
            model Outcome_obs = trt VISIT trt*VISIT / ddfm=kr;
            repeated VISIT / sub = PatientID type = cs;
            /* 周辺平均と群間差の推定 */
            estimate 'Week 36 Effect' trt 1 -1 trt*VISIT 0 0 1 0 0 -1 / cl alpha=0.05 e;
            estimate 'Week 24 Effect' trt 1 -1 trt*VISIT 0 1 0 0 -1 0 / cl alpha=0.05 e;
            estimate 'Week 12 Effect' trt 1 -1 trt*VISIT 1 0 0 -1 0 0 / cl alpha=0.05 e;
            ods output Estimates=current_est;
        run;

        /* 4. 結果の保存 (各ループのk1, k2の情報を付与して積み上げ) */
        data current_est;
            length k1_val k2_val 8;
            set current_est(where = (Label = "Week 36 Effect"));
            k1_val = &k1;
            k2_val = &k2;
        run;

        proc append base=all_estimates data=current_est force;
        run;

    %end;
%end;

%mend tipping_mi_analysis;

%tipping_mi_analysis;

options source notes;

proc sort data=all_estimates out = all_estimates_sort;
    by k1_val k2_val;
run;

/*Rubin's ruleによる結果の統合*/
proc mianalyze data=all_estimates_sort;
    by k1_val k2_val;
    modeleffects Estimate;
    stderr StdErr;
    ods output ParameterEstimates=TP_results;
run;
/*シフトパラメータ13×13=169レコードのデータセット*/

/*有意か非有意のフラグ作成*/
data TP_results2;
    set TP_results;
    if Probt < 0.05 then Sig = 1;
    else Sig = 0;
run;

/*
今回のシフトパラメータの範囲では全て有意でした。
実務上は、有意ではなくなった点を参考に臨床的な観点などから検討を行っていくことになるかと思います。
*/
