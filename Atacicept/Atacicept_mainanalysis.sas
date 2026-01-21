/*
本ページではSAP1.0版の主要な解析の大まかな再現を行っていきます。
流れは以下になります。

中間事象および欠測への対処方針(再掲)

Step1：中間事象に対する方針
(1)治療の中止：治療方針ストラテジーに基づき、中止後のデータも含めて観測されているデータを使用
(2)レスキュー薬や療法の使用：治療方針ストラテジーに基づき、使用後のデータも含めて観測されているデータを使用
(3)併用禁止薬や療法の使用：治療方針ストラテジーに基づき、使用後のデータも含めて観測されているデータを使用
(4)透析、腎移植、または死亡：これらの中間事象後のデータは解釈が難しいため、解析から除外

Step2：間欠的(Intermittent)な欠測への対処
MARを仮定し、治療群別に多重代入(50回)
- SASのPROC MI のMCMC monotone オプションを使用

Step3：中間事象による欠測への対処
・Step2で補完後のデータにおける欠測は、仮想ストラテジーに従って、Jump to Referenceにより補完
・(4)透析、腎移植、または死亡 が生じた症例については、ワーストケースアプローチ(不利な(効果が悪い)値にランダム誤差を加えた値を補完)
→jump to Referenceで補完した値をワーストケースアプローチで置換
　＊(4)の後に他の中間事象が生じた場合は(4)の対処(ワーストケースアプローチ)を優先
　＊不利な値の定義：Step1の観察データを用いて全てのvisitおよび両治療群から得られたeGFRの1%点(ランダム誤差もStep1の観察データから推定)
・(4)が生じずに試験から脱落した場合はJump to Referenceにより補完する

Step4：統計解析
MMRMを実施し、Rubin's ruleで結果を統合

*/


/*
まず、データの読み込みと欠測状況の確認を行っていきます。
*/

*データの読み込み;
proc import datafile="C:\Users\XXX\Desktop\Journal\Atacicept\analysis_data.csv" out=imported_data dbms=csv replace;/*適宜パスを設定*/
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
    by PatientID trt;
    id Week;
    var Outcome_obs;
run;

*欠測パターンの把握;
proc mi data=wide_data nimpute=0;
    var Week0 Week12 Week24 Week36;
run;

/*
出力を見ると、欠測パターンとして想定した
①間欠的な欠測
②ある時点(中間事象後を想定)以降はずっと欠測
が無事に生成できています。
*/


/*
次に、上述のStep2：間欠的(Intermittent)な欠測への対処に対してPROC MIのmcmcオプションを用いて補完 します。
*/

proc sort data = wide_data out = wide_data_sort ;
	by trt;
run;

proc mi data=wide_data_sort seed=123 nimpute=50 out=mono;
  by trt;
  var Week0 WEEK12 WEEK24 WEEK36;
  mcmc impute=monotone nbiter=1000 niter=1000;
run;

*補完後データの確認;
proc mi data=mono nimpute=0;
	by trt;
    var Week0 Week12 Week24 Week36;
run;


*出力データをみると、両群ともに間欠的な欠測が補完できています。;
*本来はburn in期間や反復回数などが十分か確認すべきところですが、今回は省略します。;


/*
次に、先述(別ページ)のStep3：Step2で補完後のデータにおける欠測に対して、仮想ストラテジーに従って、Jump to Referenceにより補完します。
*/

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
SAP通りであればこの後に腎移植や死亡が発生した症例についてはワーストケースアプローチによる補完を行うのですが、今回の目標はEstimand関連の欠測対処法の勉強のため割愛します。
実際に解析を行う場合は、腎移植や死亡による欠測にフラグを作成し(欠測理由の変数が別途あるかもしれません)、該当データに対して最悪値として定義した値を補完することになると思います。
＊論文のSAPに記載されているSASコードでもこのステップは省略されています。
*/


/*
次に、上述のStep4：統計解析としてMMRMを実施し、Rubin's ruleで結果を統合します
*/

*MMRMを行うため縦型に変換;
proc sort data=outmi;
    by _imputation_ PatientID trt;
run;

proc transpose data=outmi out=outmi_long name=VISIT;
    by _imputation_ PatientID trt;
    var WEEK0 WEEK12 WEEK24 WEEK36;
run;

/*分散共分散構造がUNではほとんどのデータセットで収束しなかったためCSで解析を実施*/
proc mixed data=outmi_long(where=(VISIT ne "Week0")) method=reml;
	by _imputation_;
	class trt(ref = "Placebo") PatientID VISIT ;
	model Outcome_obs =  trt VISIT trt*VISIT /ddfm=kr;
	repeated VISIT/ sub = PatientID type = cs;
	lsmean trt*VISIT /cl alpha=0.05 ;
	estimate 'Week 36 Effect' trt 1 -1 trt*VISIT 0 0 1 0 0 -1 / cl alpha=0.05 e;
	estimate 'Week 24 Effect' trt 1 -1 trt*VISIT 0 1 0 0 -1 0 / cl alpha=0.05 e;
	estimate 'Week 12 Effect' trt 1 -1 trt*VISIT 1 0 0 -1 0 0 / cl alpha=0.05 e;
	ods output Estimates=Estimates;
run;


/*結果の統合*/
ods select all;
proc mianalyze data=estimates;
	where Label = "Week 36 Effect";
	modeleffects estimate;
	stderr stderr;
	ods output parameterestimates=diffs1;
run;

/*(論文と値が異なりますが)Ataciceptの有効性を示す結果が得られました*/



