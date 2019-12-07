(* Content-type: application/vnd.wolfram.cdf.text *)

(*** Wolfram CDF File ***)
(* http://www.wolfram.com/cdf *)

(* CreatedBy='Mathematica 11.2' *)

(***************************************************************************)
(*                                                                         *)
(*                                                                         *)
(*  Under the Wolfram FreeCDF terms of use, this file and its content are  *)
(*  bound by the Creative Commons BY-SA Attribution-ShareAlike license.    *)
(*                                                                         *)
(*        For additional information concerning CDF licensing, see:        *)
(*                                                                         *)
(*         www.wolfram.com/cdf/adopting-cdf/licensing-options.html         *)
(*                                                                         *)
(*                                                                         *)
(***************************************************************************)

(*CacheID: 234*)
(* Internal cache information:
NotebookFileLineBreakTest
NotebookFileLineBreakTest
NotebookDataPosition[      1088,         20]
NotebookDataLength[    157823,       3721]
NotebookOptionsPosition[    148473,       3535]
NotebookOutlinePosition[    148893,       3551]
CellTagsIndexPosition[    148850,       3548]
WindowFrame->Normal*)

(* Beginning of Notebook Content *)
Notebook[{

Cell[CellGroupData[{
Cell["\<\
Theoretical 1 D model for the impact of convection in the upper atmosphere of \
the Asian monsoon\
\>", "Title",
 CellChangeTimes->{{3.7764280567884893`*^9, 3.776428071274796*^9}, {
  3.776428183970026*^9, 
  3.776428221716158*^9}},ExpressionUUID->"4b842629-12ca-4fb3-9916-\
f4cd004afc45"],

Cell["\<\
B. Legras & S. Bucci
Laboratoire de M\[EAcute]t\[EAcute]orologie Dynamique, Paris, France\
\>", "Subtitle",
 CellChangeTimes->{{3.7805097033304577`*^9, 
  3.7805097355473967`*^9}},ExpressionUUID->"208b26fc-afe1-46d7-a92d-\
abe32610fef1"],

Cell["\<\
We use here the solution to the simplified model of advection with leakage \
with no diffusion in the vertical\
\>", "Text",
 CellChangeTimes->{{3.776428079204817*^9, 3.776428140762403*^9}, 
   3.7804203165703096`*^9},ExpressionUUID->"c9cf37f4-8fc4-41b6-bb80-\
29a53e45be7a"],

Cell["\<\
Basic solution valid for all \[Theta] and times
It is assumed that the heating rate is 0 in \[Theta]0 and grows linearly from \
\[Theta]0 to \[Theta]1 with a slope \[CapitalLambda]. It is constant equal to \
A = \[CapitalLambda] (\[Theta]1 - \[Theta]0) above.
    Assuming the source is S = s0 Exp[-\[Beta] (\[CapitalTheta] - \[Theta]0)]\
\>", "Text",
 CellChangeTimes->{{3.7764282305881*^9, 3.776428246297475*^9}, {
  3.776428711590171*^9, 3.776428850189211*^9}, {3.7765080320480995`*^9, 
  3.776508210332955*^9}},ExpressionUUID->"0064d11b-757f-4896-93f6-\
db41f8147880"],

Cell[CellGroupData[{

Cell["Basic formula of the impact response", "Chapter",
 CellChangeTimes->{{3.7765082190000086`*^9, 
  3.7765082563204927`*^9}},ExpressionUUID->"00369a74-86ec-41e3-88b7-\
4a6f06ad123d"],

Cell[TextData[{
 "If \[Theta] < \[Theta]1 : F[\[Theta], t] = Exp[-(\[Alpha] + \
\[CapitalLambda]) t - \[Beta] (\[Theta] - \[Theta]0) Exp[-\[CapitalLambda] \
t]]\n\nIf \[Theta]>\[Theta]1 and ",
 Cell[BoxData[
  RowBox[{"t", "<", 
   FractionBox[
    RowBox[{"\[Theta]", "-", "\[Theta]1"}], "A"]}]],
  InitializationCell->True,
  CellChangeTimes->{{3.776405537850178*^9, 3.7764057617064676`*^9}, {
   3.77640608030727*^9, 3.776406082955756*^9}, {3.776406136981933*^9, 
   3.7764061469082727`*^9}, {3.77640684669558*^9, 3.7764068546770287`*^9}, {
   3.776406958311306*^9, 3.776406971223948*^9}, {3.776407277964577*^9, 
   3.7764072981730504`*^9}, {3.776428276354994*^9, 3.776428287832181*^9}, {
   3.776432852850367*^9, 3.776432909943634*^9}, {3.776435332485351*^9, 
   3.77643567009094*^9}, {3.776436966185503*^9, 3.776436967889153*^9}, {
   3.7764370223209543`*^9, 3.7764370486559057`*^9}, {3.7764375370627737`*^9, 
   3.776437549404471*^9}},ExpressionUUID->
  "4e6f3716-ea94-4de0-bccf-1bce04ed0899"],
 ": F[\[Theta], t] = Exp[-(\[Alpha] - \[Beta] A)t \
-\[Beta](\[Theta]-\[Theta]0)]\n\nIf \[Theta]>\[Theta]1 and ",
 Cell[BoxData[
  RowBox[{"t", ">", 
   FractionBox[
    RowBox[{"\[Theta]", "-", "\[Theta]1"}], "A"]}]],
  InitializationCell->True,
  CellChangeTimes->{{3.776405537850178*^9, 3.7764057617064676`*^9}, {
   3.77640608030727*^9, 3.776406082955756*^9}, {3.776406136981933*^9, 
   3.7764061469082727`*^9}, {3.77640684669558*^9, 3.7764068546770287`*^9}, {
   3.776406958311306*^9, 3.776406971223948*^9}, {3.776407277964577*^9, 
   3.7764072981730504`*^9}, {3.776428276354994*^9, 3.776428287832181*^9}, {
   3.776432852850367*^9, 3.776432909943634*^9}, {3.776435332485351*^9, 
   3.77643567009094*^9}, {3.776436966185503*^9, 3.776436967889153*^9}, {
   3.7764370223209543`*^9, 3.7764370486559057`*^9}, {3.7764375370627737`*^9, 
   3.776437549404471*^9}},ExpressionUUID->
  "f939f79b-1a76-40a6-a7d8-4900a862fb8a"],
 ":  F[\[Theta], t] = ",
 Cell[BoxData[
  RowBox[{"Exp", "[", 
   RowBox[{
    RowBox[{
     RowBox[{"-", "\[Alpha]"}], 
     FractionBox[
      RowBox[{"\[Theta]", "-", "\[Theta]1"}], "A"]}], "-", 
    RowBox[{
     RowBox[{"(", 
      RowBox[{"\[Alpha]", " ", "+", "\[CapitalLambda]"}], ")"}], 
     RowBox[{"(", 
      RowBox[{"t", "-", 
       FractionBox[
        RowBox[{"\[Theta]", "-", "\[Theta]1"}], "A"]}], ")"}]}], " ", "-", 
    RowBox[{
     FractionBox[
      RowBox[{"\[Beta]", " ", "A"}], "\[CapitalLambda]"], 
     RowBox[{"Exp", "[", 
      RowBox[{
       RowBox[{"-", "\[CapitalLambda]"}], 
       RowBox[{"(", 
        RowBox[{"t", "-", 
         FractionBox[
          RowBox[{"\[Theta]", "-", "\[Theta]1"}], "A"]}], ")"}]}], "]"}]}]}], 
   "]"}]],
  InitializationCell->True,
  CellChangeTimes->{{3.776405537850178*^9, 3.7764057617064676`*^9}, {
   3.77640608030727*^9, 3.776406082955756*^9}, {3.776406136981933*^9, 
   3.7764061469082727`*^9}, {3.77640684669558*^9, 3.7764068546770287`*^9}, {
   3.776406958311306*^9, 3.776406971223948*^9}, {3.776407277964577*^9, 
   3.7764072981730504`*^9}, {3.776428276354994*^9, 3.776428287832181*^9}, {
   3.776432852850367*^9, 3.776432909943634*^9}, {3.776435332485351*^9, 
   3.77643567009094*^9}, {3.776436966185503*^9, 3.776436967889153*^9}, {
   3.7764370223209543`*^9, 3.7764370486559057`*^9}, {3.7764375370627737`*^9, 
   3.776437549404471*^9}},ExpressionUUID->
  "321bc23d-cce9-4845-9a5c-40feba946e69"],
 "=\n"
}], "Text",
 CellChangeTimes->{{3.7803758318722553`*^9, 3.7803760522405987`*^9}, {
  3.780479266522069*^9, 
  3.780479334179785*^9}},ExpressionUUID->"7781d49e-ea5e-4269-b06f-\
d6bef7d370d1"],

Cell[BoxData[
 RowBox[{
  RowBox[{"F", "[", 
   RowBox[{"\[Theta]_", ",", "t_"}], "]"}], " ", ":=", " ", 
  RowBox[{"If", "[", 
   RowBox[{
    RowBox[{"\[Theta]", ">", "\[Theta]1"}], ",", 
    RowBox[{"If", "[", 
     RowBox[{
      RowBox[{"t", ">", 
       FractionBox[
        RowBox[{"\[Theta]", "-", "\[Theta]1"}], "A"]}], ",", 
      RowBox[{"Exp", "[", 
       RowBox[{
        RowBox[{
         RowBox[{"-", "\[Alpha]"}], " ", "t"}], "-", 
        RowBox[{"\[CapitalLambda]", 
         RowBox[{"(", 
          RowBox[{"t", "-", 
           FractionBox[
            RowBox[{"\[Theta]", "-", "\[Theta]1"}], "A"]}], ")"}]}], " ", "-", 
        RowBox[{
         FractionBox[
          RowBox[{"\[Beta]", " ", "A"}], "\[CapitalLambda]"], 
         RowBox[{"Exp", "[", 
          RowBox[{
           RowBox[{"-", "\[CapitalLambda]"}], 
           RowBox[{"(", 
            RowBox[{"t", "-", 
             FractionBox[
              RowBox[{"\[Theta]", "-", "\[Theta]1"}], "A"]}], ")"}]}], 
          "]"}]}]}], "]"}], ",", 
      RowBox[{"Exp", "[", 
       RowBox[{
        RowBox[{
         RowBox[{"-", 
          RowBox[{"(", 
           RowBox[{"\[Alpha]", " ", "-", " ", 
            RowBox[{"\[Beta]", " ", "A"}]}], ")"}]}], "t"}], " ", "-", 
        RowBox[{"\[Beta]", 
         RowBox[{"(", 
          RowBox[{"\[Theta]", "-", "\[Theta]0"}], ")"}]}]}], "]"}]}], "]"}], 
    ",", 
    RowBox[{"Exp", "[", 
     RowBox[{
      RowBox[{
       RowBox[{"-", 
        RowBox[{"(", 
         RowBox[{"\[Alpha]", "+", "\[CapitalLambda]"}], ")"}]}], "t"}], " ", 
      "-", " ", 
      RowBox[{"\[Beta]", 
       RowBox[{"(", 
        RowBox[{"\[Theta]", "-", "\[Theta]0"}], ")"}], " ", 
       RowBox[{"Exp", "[", 
        RowBox[{
         RowBox[{"-", "\[CapitalLambda]"}], " ", "t"}], "]"}]}]}], "]"}]}], 
   "]"}]}]], "Input",
 InitializationCell->True,
 CellChangeTimes->{{3.776405537850178*^9, 3.7764057617064676`*^9}, {
  3.77640608030727*^9, 3.776406082955756*^9}, {3.776406136981933*^9, 
  3.7764061469082727`*^9}, {3.77640684669558*^9, 3.7764068546770287`*^9}, {
  3.776406958311306*^9, 3.776406971223948*^9}, {3.776407277964577*^9, 
  3.7764072981730504`*^9}, {3.776428276354994*^9, 3.776428287832181*^9}, {
  3.776432852850367*^9, 3.776432909943634*^9}, {3.776435332485351*^9, 
  3.77643567009094*^9}, {3.776436966185503*^9, 3.776436967889153*^9}, {
  3.7764370223209543`*^9, 3.7764370486559057`*^9}, {3.7764375370627737`*^9, 
  3.776437549404471*^9}},ExpressionUUID->"3eeefdd3-c614-4e2f-a9bc-\
ffaff7edafbb"]
}, Open  ]],

Cell[CellGroupData[{

Cell["Sections", "Chapter",
 CellChangeTimes->{{3.7765167804101777`*^9, 
  3.7765167908181953`*^9}},ExpressionUUID->"80f2d523-6d0b-4d7b-a5d6-\
971d478aafa5"],

Cell[CellGroupData[{

Cell["Section in  \[Theta] for a given time (age)", "Section",
 CellChangeTimes->{{3.776508263557572*^9, 3.776508290312253*^9}, {
  3.7765167090610356`*^9, 
  3.7765167198448753`*^9}},ExpressionUUID->"350606bb-2e08-499d-94fa-\
d090cfbd0cb0"],

Cell[CellGroupData[{

Cell[BoxData[
 RowBox[{"Manipulate", "[", 
  RowBox[{
   RowBox[{"Plot", "[", 
    RowBox[{
     RowBox[{
      RowBox[{"F", "[", 
       RowBox[{"\[CapitalTheta]", ",", "t"}], "]"}], "//.", 
      "\[VeryThinSpace]", 
      RowBox[{"{", "\[VeryThinSpace]", 
       RowBox[{
        RowBox[{"\[Theta]1", "\[Rule]", 
         RowBox[{"\[Theta]0", "+", " ", 
          FractionBox["A", "\[CapitalLambda]"]}]}], ",", 
        RowBox[{"\[Alpha]", "\[Rule]", 
         FractionBox["1", "10"]}], ",", 
        RowBox[{"\[Beta]", "\[Rule]", "b"}], ",", 
        RowBox[{"A", "\[Rule]", "a"}], ",", 
        RowBox[{"\[CapitalLambda]", "\[Rule]", "\[Lambda]"}], ",", 
        RowBox[{"\[Theta]0", "\[Rule]", "\[CurlyRho]"}]}], "}"}]}], ",", 
     RowBox[{"{", 
      RowBox[{"\[CapitalTheta]", ",", "\[CurlyRho]", ",", "420"}], "}"}], ",", 
     RowBox[{"PlotRange", "\[Rule]", "All"}], ",", 
     RowBox[{"Filling", "\[Rule]", "Bottom"}]}], "]"}], ",", 
   RowBox[{"{", 
    RowBox[{
     RowBox[{"{", 
      RowBox[{"a", ",", "1", ",", "\"\<A\>\""}], "}"}], ",", "1", ",", "1.5", 
     ",", 
     RowBox[{"Appearance", "\[Rule]", "\"\<Labeled\>\""}]}], "}"}], ",", 
   RowBox[{"{", 
    RowBox[{
     RowBox[{"{", 
      RowBox[{"b", ",", "0.365", ",", "\"\<\[Beta]\>\""}], "}"}], ",", "0.1", 
     ",", "0.6", ",", 
     RowBox[{"Appearance", "\[Rule]", "\"\<Labeled\>\""}]}], "}"}], ",", 
   RowBox[{"{", 
    RowBox[{
     RowBox[{"{", 
      RowBox[{"\[Lambda]", ",", "0.1", ",", "\"\<\[CapitalLambda]\>\""}], 
      "}"}], ",", "0.1", ",", "0.2", ",", 
     RowBox[{"Appearance", "\[Rule]", "\"\<Labeled\>\""}]}], "}"}], ",", 
   RowBox[{"{", 
    RowBox[{
     RowBox[{"{", 
      RowBox[{"\[CurlyRho]", ",", "360", ",", "\"\<\[Theta]0\>\""}], "}"}], 
     ",", "350", ",", "370", ",", 
     RowBox[{"Appearance", "\[Rule]", "\"\<Labeled\>\""}]}], "}"}], ",", 
   RowBox[{"{", 
    RowBox[{
     RowBox[{"{", 
      RowBox[{"t", ",", "30", ",", "\"\<t\>\""}], "}"}], ",", "0", ",", "60", 
     ",", 
     RowBox[{"Appearance", "\[Rule]", "\"\<Labeled\>\""}]}], "}"}], ",", 
   RowBox[{"SaveDefinitions", "\[Rule]", "True"}]}], "]"}]], "Input",
 InitializationCell->True,
 CellChangeTimes->{{3.7764089115682364`*^9, 3.7764090226098704`*^9}, {
   3.776409067331358*^9, 3.7764090779762816`*^9}, {3.7764093224678173`*^9, 
   3.776409572396212*^9}, {3.7764096224396653`*^9, 3.7764096638050823`*^9}, {
   3.7764097399995947`*^9, 3.7764097578963833`*^9}, {3.776409799262787*^9, 
   3.776409833496255*^9}, {3.7764099851317463`*^9, 3.77641000493445*^9}, {
   3.7764100382314663`*^9, 3.7764100440411367`*^9}, {3.7764101549478645`*^9, 
   3.7764101610623755`*^9}, {3.7764104033728766`*^9, 3.776410436862774*^9}, {
   3.7764104795384417`*^9, 3.776410549721387*^9}, 3.776428331968004*^9, {
   3.7764371243388977`*^9, 3.7764371303715067`*^9}, 3.776437164978725*^9, {
   3.77643721536313*^9, 3.7764372240422792`*^9}, {3.776439075019648*^9, 
   3.776439140184679*^9}, {3.776439208387807*^9, 3.776439214867412*^9}, {
   3.7803093270512967`*^9, 3.7803093277856474`*^9}, 3.7803165264405375`*^9, {
   3.780340543873149*^9, 3.780340550628499*^9}, {3.7803449188859615`*^9, 
   3.7803449200109167`*^9}, {3.7805094500665026`*^9, 
   3.7805094556565866`*^9}, {3.78050949684064*^9, 3.780509514420972*^9}, {
   3.780559252258082*^9, 
   3.780559258649933*^9}},ExpressionUUID->"d5610d9c-1d00-4816-8a21-\
020cc0bc0608"],

Cell[BoxData[
 TagBox[
  StyleBox[
   DynamicModuleBox[{$CellContext`a$$ = 1, $CellContext`b$$ = 
    0.365, $CellContext`t$$ = 30, $CellContext`\[Lambda]$$ = 
    0.1, $CellContext`\[CurlyRho]$$ = 360, Typeset`show$$ = True, 
    Typeset`bookmarkList$$ = {}, Typeset`bookmarkMode$$ = "Menu", 
    Typeset`animator$$, Typeset`animvar$$ = 1, Typeset`name$$ = 
    "\"untitled\"", Typeset`specs$$ = {{{
       Hold[$CellContext`a$$], 1, "A"}, 1, 1.5}, {{
       Hold[$CellContext`b$$], 0.365, "\[Beta]"}, 0.1, 0.6}, {{
       Hold[$CellContext`\[Lambda]$$], 0.1, "\[CapitalLambda]"}, 0.1, 0.2}, {{
       Hold[$CellContext`\[CurlyRho]$$], 360, "\[Theta]0"}, 350, 370}, {{
       Hold[$CellContext`t$$], 30, "t"}, 0, 60}}, Typeset`size$$ = {
    360., {104., 109.}}, Typeset`update$$ = 0, Typeset`initDone$$, 
    Typeset`skipInitDone$$ = False, $CellContext`a$18319$$ = 
    0, $CellContext`b$18320$$ = 0, $CellContext`\[Lambda]$18321$$ = 
    0, $CellContext`\[CurlyRho]$18322$$ = 0, $CellContext`t$18323$$ = 0}, 
    DynamicBox[Manipulate`ManipulateBoxes[
     1, StandardForm, 
      "Variables" :> {$CellContext`a$$ = 1, $CellContext`b$$ = 
        0.365, $CellContext`t$$ = 30, $CellContext`\[Lambda]$$ = 
        0.1, $CellContext`\[CurlyRho]$$ = 360}, "ControllerVariables" :> {
        Hold[$CellContext`a$$, $CellContext`a$18319$$, 0], 
        Hold[$CellContext`b$$, $CellContext`b$18320$$, 0], 
        Hold[$CellContext`\[Lambda]$$, $CellContext`\[Lambda]$18321$$, 0], 
        Hold[$CellContext`\[CurlyRho]$$, $CellContext`\[CurlyRho]$18322$$, 0], 
        Hold[$CellContext`t$$, $CellContext`t$18323$$, 0]}, 
      "OtherVariables" :> {
       Typeset`show$$, Typeset`bookmarkList$$, Typeset`bookmarkMode$$, 
        Typeset`animator$$, Typeset`animvar$$, Typeset`name$$, 
        Typeset`specs$$, Typeset`size$$, Typeset`update$$, Typeset`initDone$$,
         Typeset`skipInitDone$$}, "Body" :> Plot[
        ReplaceRepeated[
         $CellContext`F[$CellContext`\[CapitalTheta], $CellContext`t$$], \
{$CellContext`\[Theta]1 -> $CellContext`\[Theta]0 + \
$CellContext`A/$CellContext`\[CapitalLambda], $CellContext`\[Alpha] -> 
          1/10, $CellContext`\[Beta] -> $CellContext`b$$, $CellContext`A -> \
$CellContext`a$$, $CellContext`\[CapitalLambda] -> $CellContext`\[Lambda]$$, \
$CellContext`\[Theta]0 -> $CellContext`\[CurlyRho]$$}], {$CellContext`\
\[CapitalTheta], $CellContext`\[CurlyRho]$$, 420}, PlotRange -> All, Filling -> 
        Bottom], 
      "Specifications" :> {{{$CellContext`a$$, 1, "A"}, 1, 1.5, Appearance -> 
         "Labeled"}, {{$CellContext`b$$, 0.365, "\[Beta]"}, 0.1, 0.6, 
         Appearance -> 
         "Labeled"}, {{$CellContext`\[Lambda]$$, 0.1, "\[CapitalLambda]"}, 
         0.1, 0.2, Appearance -> 
         "Labeled"}, {{$CellContext`\[CurlyRho]$$, 360, "\[Theta]0"}, 350, 
         370, Appearance -> "Labeled"}, {{$CellContext`t$$, 30, "t"}, 0, 60, 
         Appearance -> "Labeled"}}, "Options" :> {}, "DefaultOptions" :> {}],
     ImageSizeCache->{411., {207., 213.}},
     SingleEvaluation->True],
    Deinitialization:>None,
    DynamicModuleValues:>{},
    Initialization:>({$CellContext`F[
         Pattern[$CellContext`\[Theta], 
          Blank[]], 
         Pattern[$CellContext`t, 
          Blank[]]] := If[$CellContext`\[Theta] > $CellContext`\[Theta]1, 
         If[$CellContext`t > ($CellContext`\[Theta] - \
$CellContext`\[Theta]1)/$CellContext`A, 
          
          Exp[(-$CellContext`\[Alpha]) $CellContext`t - $CellContext`\
\[CapitalLambda] ($CellContext`t - ($CellContext`\[Theta] - $CellContext`\
\[Theta]1)/$CellContext`A) - (($CellContext`\[Beta] \
$CellContext`A)/$CellContext`\[CapitalLambda]) 
           Exp[(-$CellContext`\[CapitalLambda]) ($CellContext`t - \
($CellContext`\[Theta] - $CellContext`\[Theta]1)/$CellContext`A)]], 
          
          Exp[(-($CellContext`\[Alpha] - $CellContext`\[Beta] \
$CellContext`A)) $CellContext`t - $CellContext`\[Beta] ($CellContext`\[Theta] - \
$CellContext`\[Theta]0)]], 
         Exp[(-($CellContext`\[Alpha] + $CellContext`\[CapitalLambda])) \
$CellContext`t - ($CellContext`\[Beta] ($CellContext`\[Theta] - $CellContext`\
\[Theta]0)) Exp[(-$CellContext`\[CapitalLambda]) $CellContext`t]]]}; 
     Typeset`initDone$$ = True),
    SynchronousInitialization->True,
    UndoTrackedVariables:>{Typeset`show$$, Typeset`bookmarkMode$$},
    UnsavedVariables:>{Typeset`initDone$$},
    UntrackedVariables:>{Typeset`size$$}], "Manipulate",
   Deployed->True,
   StripOnInput->False],
  Manipulate`InterpretManipulate[1]]], "Output",
 CellChangeTimes->{
  3.776455575393448*^9, 3.7765080258276634`*^9, 3.7765812178014793`*^9, 
   3.7767066317122173`*^9, 3.776707174589498*^9, {3.776707284309761*^9, 
   3.776707335253644*^9}, 3.7767074005930014`*^9, 3.7803091651132736`*^9, 
   3.780309335628809*^9, 3.7803174786280327`*^9, 3.7803405549631777`*^9, 
   3.7803449299033566`*^9, {3.7804223087952633`*^9, 3.780422338432567*^9}, 
   3.7805094573908587`*^9, 3.7805095168582954`*^9, {3.7805144279189167`*^9, 
   3.7805144331529226`*^9}, 3.7805149276344414`*^9, 
   3.780560759875751*^9},ExpressionUUID->"70e2c71a-58e6-448b-ba15-\
5b0190dd85ea"]
}, Open  ]]
}, Open  ]],

Cell[CellGroupData[{

Cell["Section in t (age) for a given \[Theta]", "Section",
 CellChangeTimes->{{3.776508263557572*^9, 
  3.776508290312253*^9}},ExpressionUUID->"6b94b792-a36b-4cd3-bc80-\
6ad9d9535a24"],

Cell[CellGroupData[{

Cell[BoxData[
 RowBox[{"Manipulate", "[", 
  RowBox[{
   RowBox[{"Plot", "[", 
    RowBox[{
     RowBox[{
      RowBox[{"F", "[", 
       RowBox[{"\[CapitalTheta]", ",", "t"}], "]"}], "//.", 
      "\[VeryThinSpace]", 
      RowBox[{"{", "\[VeryThinSpace]", 
       RowBox[{
        RowBox[{"\[Theta]1", "\[Rule]", 
         RowBox[{"\[Theta]0", "+", " ", 
          FractionBox["A", "\[CapitalLambda]"]}]}], ",", 
        RowBox[{"\[Alpha]", "\[Rule]", 
         FractionBox["1", "10"]}], ",", 
        RowBox[{"\[Beta]", "\[Rule]", "b"}], ",", 
        RowBox[{"A", "\[Rule]", "a"}], ",", 
        RowBox[{"\[CapitalLambda]", "\[Rule]", "\[Lambda]"}], ",", 
        RowBox[{"\[Theta]0", "\[Rule]", "\[CurlyRho]"}]}], "}"}]}], ",", 
     RowBox[{"{", 
      RowBox[{"t", ",", "0", ",", "60"}], "}"}], ",", 
     RowBox[{"PlotRange", "\[Rule]", "All"}], ",", 
     RowBox[{"Filling", "\[Rule]", "Bottom"}]}], "]"}], ",", 
   RowBox[{"{", 
    RowBox[{
     RowBox[{"{", 
      RowBox[{"a", ",", "1", ",", "\"\<A\>\""}], "}"}], ",", "1", ",", "1.5", 
     ",", 
     RowBox[{"Appearance", "\[Rule]", "\"\<Labeled\>\""}]}], "}"}], ",", 
   RowBox[{"{", 
    RowBox[{
     RowBox[{"{", 
      RowBox[{"b", ",", "0.365", ",", "\"\<\[Beta]\>\""}], "}"}], ",", "0.05",
      ",", "0.6", ",", 
     RowBox[{"Appearance", "\[Rule]", "\"\<Labeled\>\""}]}], "}"}], ",", 
   RowBox[{"{", 
    RowBox[{
     RowBox[{"{", 
      RowBox[{"\[Lambda]", ",", "0.1", ",", "\"\<\[CapitalLambda]\>\""}], 
      "}"}], ",", "0.1", ",", "0.2", ",", 
     RowBox[{"Appearance", "\[Rule]", "\"\<Labeled\>\""}]}], "}"}], ",", 
   RowBox[{"{", 
    RowBox[{
     RowBox[{"{", 
      RowBox[{"\[CurlyRho]", ",", "360", ",", "\"\<\[Theta]0\>\""}], "}"}], 
     ",", "350", ",", "370", ",", 
     RowBox[{"Appearance", "\[Rule]", "\"\<Labeled\>\""}]}], "}"}], ",", 
   RowBox[{"{", 
    RowBox[{
     RowBox[{"{", 
      RowBox[{"\[CapitalTheta]", ",", "375", ",", "\"\<\[Theta]\>\""}], "}"}],
      ",", "\[CurlyRho]", ",", "420", ",", 
     RowBox[{"Appearance", "\[Rule]", "\"\<Labeled\>\""}]}], "}"}], ",", 
   RowBox[{"SaveDefinitions", "\[Rule]", "True"}]}], "]"}]], "Input",
 InitializationCell->True,
 CellChangeTimes->{{3.7764089115682364`*^9, 3.7764090226098704`*^9}, {
   3.776409067331358*^9, 3.7764090779762816`*^9}, {3.7764093224678173`*^9, 
   3.776409572396212*^9}, {3.7764096224396653`*^9, 3.7764096638050823`*^9}, {
   3.7764097399995947`*^9, 3.7764097578963833`*^9}, {3.776409799262787*^9, 
   3.776409833496255*^9}, {3.7764099851317463`*^9, 3.77641000493445*^9}, {
   3.7764100382314663`*^9, 3.7764100440411367`*^9}, {3.7764101549478645`*^9, 
   3.7764101610623755`*^9}, {3.7764104033728766`*^9, 3.776410436862774*^9}, {
   3.7764104795384417`*^9, 3.776410549721387*^9}, 3.776428331968004*^9, {
   3.7764371243388977`*^9, 3.7764371303715067`*^9}, 3.776437164978725*^9, {
   3.77643721536313*^9, 3.7764372240422792`*^9}, 3.7767099180717344`*^9, {
   3.776711402405743*^9, 3.7767114043386393`*^9}, {3.776711435377899*^9, 
   3.7767114392007113`*^9}, {3.780316230390089*^9, 3.780316230944525*^9}, {
   3.78031650240337*^9, 3.7803165328993626`*^9}, {3.780340567151187*^9, 
   3.7803405724646635`*^9}, {3.7803425400397644`*^9, 3.780342541850659*^9}, {
   3.7803448252763586`*^9, 3.780344826953353*^9}, {3.780509553819659*^9, 
   3.7805095953137636`*^9}, {3.7805592680777726`*^9, 
   3.780559274907071*^9}},ExpressionUUID->"275c7960-e4da-44b9-ba32-\
b10e12500197"],

Cell[BoxData[
 TagBox[
  StyleBox[
   DynamicModuleBox[{$CellContext`a$$ = 1, $CellContext`b$$ = 
    0.365, $CellContext`\[CapitalTheta]$$ = 375, $CellContext`\[Lambda]$$ = 
    0.1, $CellContext`\[CurlyRho]$$ = 360, Typeset`show$$ = True, 
    Typeset`bookmarkList$$ = {}, Typeset`bookmarkMode$$ = "Menu", 
    Typeset`animator$$, Typeset`animvar$$ = 1, Typeset`name$$ = 
    "\"untitled\"", Typeset`specs$$ = {{{
       Hold[$CellContext`a$$], 1, "A"}, 1, 1.5}, {{
       Hold[$CellContext`b$$], 0.365, "\[Beta]"}, 0.05, 0.6}, {{
       Hold[$CellContext`\[Lambda]$$], 0.1, "\[CapitalLambda]"}, 0.1, 0.2}, {{
       Hold[$CellContext`\[CurlyRho]$$], 360, "\[Theta]0"}, 350, 370}, {{
       Hold[$CellContext`\[CapitalTheta]$$], 375, "\[Theta]"}, 
      Dynamic[$CellContext`\[CurlyRho]$$], 420}}, Typeset`size$$ = {
    360., {105., 109.}}, Typeset`update$$ = 0, Typeset`initDone$$, 
    Typeset`skipInitDone$$ = False, $CellContext`a$18406$$ = 
    0, $CellContext`b$18407$$ = 0, $CellContext`\[Lambda]$18408$$ = 
    0, $CellContext`\[CurlyRho]$18409$$ = 
    0, $CellContext`\[CapitalTheta]$18410$$ = 0}, 
    DynamicBox[Manipulate`ManipulateBoxes[
     1, StandardForm, 
      "Variables" :> {$CellContext`a$$ = 1, $CellContext`b$$ = 
        0.365, $CellContext`\[CapitalTheta]$$ = 375, $CellContext`\[Lambda]$$ = 
        0.1, $CellContext`\[CurlyRho]$$ = 360}, "ControllerVariables" :> {
        Hold[$CellContext`a$$, $CellContext`a$18406$$, 0], 
        Hold[$CellContext`b$$, $CellContext`b$18407$$, 0], 
        Hold[$CellContext`\[Lambda]$$, $CellContext`\[Lambda]$18408$$, 0], 
        Hold[$CellContext`\[CurlyRho]$$, $CellContext`\[CurlyRho]$18409$$, 0], 
        Hold[$CellContext`\[CapitalTheta]$$, \
$CellContext`\[CapitalTheta]$18410$$, 0]}, 
      "OtherVariables" :> {
       Typeset`show$$, Typeset`bookmarkList$$, Typeset`bookmarkMode$$, 
        Typeset`animator$$, Typeset`animvar$$, Typeset`name$$, 
        Typeset`specs$$, Typeset`size$$, Typeset`update$$, Typeset`initDone$$,
         Typeset`skipInitDone$$}, "Body" :> Plot[
        ReplaceRepeated[
         $CellContext`F[$CellContext`\[CapitalTheta]$$, $CellContext`t], \
{$CellContext`\[Theta]1 -> $CellContext`\[Theta]0 + \
$CellContext`A/$CellContext`\[CapitalLambda], $CellContext`\[Alpha] -> 
          1/10, $CellContext`\[Beta] -> $CellContext`b$$, $CellContext`A -> \
$CellContext`a$$, $CellContext`\[CapitalLambda] -> $CellContext`\[Lambda]$$, \
$CellContext`\[Theta]0 -> $CellContext`\[CurlyRho]$$}], {$CellContext`t, 0, 
         60}, PlotRange -> All, Filling -> Bottom], 
      "Specifications" :> {{{$CellContext`a$$, 1, "A"}, 1, 1.5, Appearance -> 
         "Labeled"}, {{$CellContext`b$$, 0.365, "\[Beta]"}, 0.05, 0.6, 
         Appearance -> 
         "Labeled"}, {{$CellContext`\[Lambda]$$, 0.1, "\[CapitalLambda]"}, 
         0.1, 0.2, Appearance -> 
         "Labeled"}, {{$CellContext`\[CurlyRho]$$, 360, "\[Theta]0"}, 350, 
         370, Appearance -> 
         "Labeled"}, {{$CellContext`\[CapitalTheta]$$, 375, "\[Theta]"}, 
         Dynamic[$CellContext`\[CurlyRho]$$], 420, Appearance -> "Labeled"}}, 
      "Options" :> {}, "DefaultOptions" :> {}],
     ImageSizeCache->{411., {207., 213.}},
     SingleEvaluation->True],
    Deinitialization:>None,
    DynamicModuleValues:>{},
    Initialization:>({$CellContext`F[
         Pattern[$CellContext`\[Theta], 
          Blank[]], 
         Pattern[$CellContext`t, 
          Blank[]]] := If[$CellContext`\[Theta] > $CellContext`\[Theta]1, 
         If[$CellContext`t > ($CellContext`\[Theta] - \
$CellContext`\[Theta]1)/$CellContext`A, 
          
          Exp[(-$CellContext`\[Alpha]) $CellContext`t - $CellContext`\
\[CapitalLambda] ($CellContext`t - ($CellContext`\[Theta] - $CellContext`\
\[Theta]1)/$CellContext`A) - (($CellContext`\[Beta] \
$CellContext`A)/$CellContext`\[CapitalLambda]) 
           Exp[(-$CellContext`\[CapitalLambda]) ($CellContext`t - \
($CellContext`\[Theta] - $CellContext`\[Theta]1)/$CellContext`A)]], 
          
          Exp[(-($CellContext`\[Alpha] - $CellContext`\[Beta] \
$CellContext`A)) $CellContext`t - $CellContext`\[Beta] ($CellContext`\[Theta] - \
$CellContext`\[Theta]0)]], 
         Exp[(-($CellContext`\[Alpha] + $CellContext`\[CapitalLambda])) \
$CellContext`t - ($CellContext`\[Beta] ($CellContext`\[Theta] - $CellContext`\
\[Theta]0)) Exp[(-$CellContext`\[CapitalLambda]) $CellContext`t]]]}; 
     Typeset`initDone$$ = True),
    SynchronousInitialization->True,
    UndoTrackedVariables:>{Typeset`show$$, Typeset`bookmarkMode$$},
    UnsavedVariables:>{Typeset`initDone$$},
    UntrackedVariables:>{Typeset`size$$}], "Manipulate",
   Deployed->True,
   StripOnInput->False],
  Manipulate`InterpretManipulate[1]]], "Output",
 CellChangeTimes->{3.780422309014*^9, 3.7804223798367815`*^9, 
  3.7805095659338236`*^9, 3.7805095969386115`*^9, 3.7805144281532784`*^9, 
  3.780514927931293*^9, 
  3.7805607602975955`*^9},ExpressionUUID->"2c66fc81-5a79-4299-9919-\
96282d97e0f6"]
}, Open  ]],

Cell[TextData[{
 "For \[Theta] < \[Theta]1, the derivative in t = 0 is proportional and of \
the sign  of - (\[Alpha] + \[CapitalLambda]) + \[Beta] (\[Theta] - \
\[Theta]0). It is positive for  \[Theta] - \[Theta]0 > ",
 Cell[BoxData[
  FormBox[
   FractionBox[
    RowBox[{"\[Alpha]", "+", "\[CapitalLambda]"}], 
    RowBox[{"\[Beta]", " ", "\[CapitalLambda]"}]], TraditionalForm]],
  ExpressionUUID->"826a74fb-9cc8-4481-a5e5-3010cfaf005e"],
 "which occurs only if  ",
 Cell[BoxData[
  FormBox[
   FractionBox[
    RowBox[{"\[Alpha]", "+", "\[CapitalLambda]"}], 
    RowBox[{"\[Beta]", " ", "\[CapitalLambda]"}]], TraditionalForm]],
  ExpressionUUID->"5eec2015-459f-4d6a-a268-eaca2dc5d5b0"],
 " < ",
 Cell[BoxData[
  FormBox[
   FractionBox["A", "\[CapitalLambda]"], TraditionalForm]],ExpressionUUID->
  "2785a2cb-b907-4a1f-a398-582099f37a22"],
 ", that is \[Beta] > ",
 Cell[BoxData[
  FormBox[
   FractionBox[
    RowBox[{"\[Alpha]", "+", "\[CapitalLambda]"}], "A"], TraditionalForm]],
  ExpressionUUID->"7a1402d4-5caf-4a78-968f-6a23b9a047ba"],
 ". It is otherwise negative. When the derivative in t=0 is positive, there \
is a maximum as a function of time which is obtained  for t = ",
 Cell[BoxData[
  FormBox[
   RowBox[{
    FractionBox["1", "\[CapitalLambda]"], 
    RowBox[{"Log", "(", 
     FractionBox[
      RowBox[{"\[Beta]", " ", "\[CapitalLambda]", " ", 
       RowBox[{"(", 
        RowBox[{"\[Theta]", "-", "\[Theta]0"}], ")"}]}], 
      RowBox[{"\[Alpha]", "+", "\[CapitalLambda]"}]]}]}], TraditionalForm]],
  ExpressionUUID->"b3a8a3d0-a767-44e2-877f-d96b6b82dc69"],
 ")\nFor \[Theta] > \[Theta]1, the derivative in t=0 is proportional to \
\[Beta] A - \[Alpha], it is positive if \[Beta]> \[Alpha]/A. If \[Beta] = \
\[Alpha]/A,  F does not depend on t over 0 < t < (\[Theta] - \[Theta]1)/A. \
For all the values of the parameters, the time derivative has same sign over \
the time interval ",
 Cell[BoxData[
  FormBox[
   RowBox[{"[", 
    RowBox[{"0", ",", 
     FractionBox[
      RowBox[{"\[Theta]", "-", "\[Theta]1"}], "A"]}], "]"}], 
   TraditionalForm]],ExpressionUUID->"0eae88ea-780b-4a40-88da-1a16b20971f4"],
 ".\nFor t = ",
 Cell[BoxData[
  FormBox[
   FractionBox[
    RowBox[{"\[Theta]", "-", "\[Theta]1"}], "A"], TraditionalForm]],
  ExpressionUUID->"86d2966e-9a60-4fce-abb3-669b19318eb6"],
 ", the right time derivative is proportional to  \
-(\[Alpha]+\[CapitalLambda])+\[Beta]A. This derivative is positive for  \
\[Beta] > ",
 Cell[BoxData[
  FormBox[
   FractionBox[
    RowBox[{"\[Alpha]", "+", "\[CapitalLambda]"}], "A"], TraditionalForm]],
  ExpressionUUID->"7e66114d-e56e-4822-b2d8-af6f4c5975c3"],
 ", leading to a smooth maximum at larger time. When \[Beta] < ",
 Cell[BoxData[
  FormBox[
   FractionBox[
    RowBox[{"\[Alpha]", "+", "\[CapitalLambda]"}], "A"], TraditionalForm]],
  ExpressionUUID->"2e5c695a-8a6e-43af-902d-4e65eac4c387"],
 ", the maximum is located at t = ",
 Cell[BoxData[
  FormBox[
   FractionBox[
    RowBox[{"\[Theta]", "-", "\[Theta]1"}], "A"], TraditionalForm]],
  ExpressionUUID->"e2dd45cd-4273-463d-aa2c-cf5bfd54092a"],
 ", and is sharp (discontinuity of the derivatives). Assuming \[Beta] > ",
 Cell[BoxData[
  FormBox[
   FractionBox[
    RowBox[{"\[Alpha]", "+", "\[CapitalLambda]"}], "A"], TraditionalForm]],
  ExpressionUUID->"d1e49e58-009d-4f45-8838-3c3ab2229b42"],
 ", the maximum occurs at t = ",
 Cell[BoxData[
  FormBox[
   FractionBox[
    RowBox[{"\[Theta]", "-", "\[Theta]1"}], "A"], TraditionalForm]],
  ExpressionUUID->"d838f0f7-f16f-422a-b042-8dce1b7abd93"],
 "+",
 Cell[BoxData[
  FormBox[
   RowBox[{
    FractionBox["1", "\[CapitalLambda]"], 
    RowBox[{"Log", "(", 
     FractionBox[
      RowBox[{"\[Beta]", " ", "A"}], 
      RowBox[{"\[Alpha]", "+", "\[CapitalLambda]"}]]}]}], TraditionalForm]],
  ExpressionUUID->"730bb580-b85f-4069-bf33-0d3445f49dd5"],
 "). Therefore, it varies linearly with \[Theta]. The large t cutoff is in \
Exp(- (\[Alpha] + \[CapitalLambda])t).  More generally, for t > ",
 Cell[BoxData[
  FormBox[
   RowBox[{
    FractionBox[
     RowBox[{"\[Theta]", "-", "\[Theta]1"}], "A"], ",", " ", 
    RowBox[{
    "adding", " ", "\[CapitalDelta]\[Theta]", " ", "to", " ", "\[Theta]", " ",
      "is", " ", "equivalent", " "}]}], TraditionalForm]],ExpressionUUID->
  "4bffc0a6-5956-429d-90cf-f0111852049d"],
 "to shifting t by \[CapitalDelta]t= ",
 Cell[BoxData[
  FormBox[
   FractionBox["\[CapitalDelta]\[Theta]", "A"], TraditionalForm]],
  ExpressionUUID->"16010bbc-ef1c-4eb5-92cd-808805753e13"],
 "and multiplying the amplitude by ",
 Cell[BoxData[
  FormBox[
   RowBox[{"Exp", "(", 
    RowBox[{"-", 
     FractionBox[
      RowBox[{"\[Alpha]", "+", "\[CapitalLambda]"}], "A"]}]}], 
   TraditionalForm]],ExpressionUUID->"4e1c796e-755e-4767-8d50-bae8b6d7feca"],
 "\[CapitalDelta]\[Theta]). For t <",
 Cell[BoxData[
  RowBox[{
   FractionBox[
    RowBox[{"\[Theta]", "-", "\[Theta]1"}], "A"], ",", " "}]],ExpressionUUID->
  "6abcb2f9-4dc6-40bc-9ade-b2c2c27b4f8e"],
 "it is equivalent to a shift \[CapitalDelta]t = ",
 Cell[BoxData[
  FormBox[
   FractionBox[
    RowBox[{"\[Beta]", " "}], 
    RowBox[{"\[Beta]A", "-", "\[Alpha]"}]], TraditionalForm]],ExpressionUUID->
  "df9d396b-b5f0-4cb1-8370-004e2cb99103"],
 "\[CapitalDelta]\[Theta] or to a multiplication of the amplitude by Exp(-\
\[Beta] \[CapitalDelta]\[Theta]).\nAs the mean age is affected by the two \
components \n"
}], "Text",
 CellChangeTimes->{{3.7803416004494696`*^9, 3.780342121675558*^9}, {
   3.7803422037727175`*^9, 3.780342206116302*^9}, {3.7803422521829515`*^9, 
   3.7803422926813593`*^9}, {3.7803423655159144`*^9, 
   3.7803423986837378`*^9}, {3.7803450314928317`*^9, 
   3.7803451574141417`*^9}, {3.7803452161436653`*^9, 
   3.7803453238148384`*^9}, {3.780345657161129*^9, 3.780345659728791*^9}, {
   3.7803458441034455`*^9, 3.78034589258246*^9}, {3.780345924400325*^9, 
   3.780345946323673*^9}, {3.7803461572062483`*^9, 3.7803462693903646`*^9}, {
   3.780346437755943*^9, 3.780346602600599*^9}, {3.780346640278222*^9, 
   3.7803466720147977`*^9}, {3.7803467080944033`*^9, 3.780346710063012*^9}, {
   3.7803467607038417`*^9, 3.7803468207433453`*^9}, {3.7803468541638875`*^9, 
   3.780346940474944*^9}, {3.780373242509693*^9, 3.7803733719758577`*^9}, {
   3.7803735092723656`*^9, 3.7803736253521824`*^9}, {3.78037371432431*^9, 
   3.780373716324165*^9}, {3.780373774593932*^9, 3.780373795002757*^9}, 
   3.7803738637619457`*^9, {3.7803739068520565`*^9, 3.7803739372538695`*^9}, {
   3.7803747484075947`*^9, 3.7803748185777664`*^9}, {3.780374990606753*^9, 
   3.780375054719812*^9}, {3.780375094910079*^9, 3.7803753075765514`*^9}, {
   3.7803753826365604`*^9, 3.7803754633169084`*^9}, {3.7803779521240935`*^9, 
   3.780377972758561*^9}, {3.78049016836825*^9, 
   3.7804901807461963`*^9}},ExpressionUUID->"464b4d9a-53aa-4de3-bc4f-\
997cb4be354b"]
}, Open  ]]
}, Open  ]],

Cell[CellGroupData[{

Cell["Age/ altitude map", "Chapter",
 CellChangeTimes->{{3.776516738299304*^9, 3.776516758663618*^9}, {
   3.780341592529708*^9, 3.7803415966841125`*^9}, 
   3.7805600993946085`*^9},ExpressionUUID->"b09a3313-ec68-4def-8a0b-\
43b2c87b1621"],

Cell[CellGroupData[{

Cell[BoxData[
 RowBox[{"Manipulate", "[", 
  RowBox[{
   RowBox[{"ContourPlot", "[", 
    RowBox[{
     RowBox[{
      RowBox[{"Log", "[", 
       RowBox[{"F", "[", 
        RowBox[{"\[CapitalTheta]", ",", "t"}], "]"}], "]"}], "//.", 
      "\[VeryThinSpace]", 
      RowBox[{"{", "\[VeryThinSpace]", 
       RowBox[{
        RowBox[{"\[Theta]1", "\[Rule]", 
         RowBox[{"\[Theta]0", "+", " ", 
          FractionBox["A", "\[CapitalLambda]"]}]}], ",", 
        RowBox[{"\[Alpha]", "\[Rule]", 
         FractionBox["1", "13"]}], ",", 
        RowBox[{"\[Beta]", "\[Rule]", "b"}], ",", 
        RowBox[{"A", "\[Rule]", "a"}], ",", 
        RowBox[{"\[CapitalLambda]", "\[Rule]", "\[Lambda]"}], ",", 
        RowBox[{"\[Theta]0", "\[Rule]", "\[CurlyRho]"}]}], "}"}]}], ",", 
     RowBox[{"{", 
      RowBox[{"t", ",", "0", ",", "60"}], "}"}], ",", 
     RowBox[{"{", 
      RowBox[{"\[CapitalTheta]", ",", "\[CurlyRho]", ",", "400"}], "}"}]}], 
    "]"}], ",", 
   RowBox[{"{", 
    RowBox[{
     RowBox[{"{", 
      RowBox[{"a", ",", "1", ",", "\"\<A\>\""}], "}"}], ",", "1", ",", "1.5", 
     ",", 
     RowBox[{"Appearance", "\[Rule]", "\"\<Labeled\>\""}]}], "}"}], ",", 
   RowBox[{"{", 
    RowBox[{
     RowBox[{"{", 
      RowBox[{"b", ",", "0.365", ",", "\"\<\[Beta]\>\""}], "}"}], ",", "0.1", 
     ",", "0.6", ",", 
     RowBox[{"Appearance", "\[Rule]", "\"\<Labeled\>\""}]}], "}"}], ",", 
   RowBox[{"{", 
    RowBox[{
     RowBox[{"{", 
      RowBox[{"\[Lambda]", ",", "0.1", ",", "\"\<\[CapitalLambda]\>\""}], 
      "}"}], ",", "0.1", ",", "0.2", ",", 
     RowBox[{"Appearance", "\[Rule]", "\"\<Labeled\>\""}]}], "}"}], ",", 
   RowBox[{"{", 
    RowBox[{
     RowBox[{"{", 
      RowBox[{"\[CurlyRho]", ",", "360", ",", "\"\<\[Theta]0\>\""}], "}"}], 
     ",", "350", ",", "370", ",", 
     RowBox[{"Appearance", "\[Rule]", "\"\<Labeled\>\""}]}], "}"}], ",", 
   RowBox[{"SaveDefinitions", "\[Rule]", "True"}]}], "]"}]], "Input",
 InitializationCell->True,
 CellChangeTimes->{{3.7764089115682364`*^9, 3.7764090226098704`*^9}, {
   3.776409067331358*^9, 3.7764090779762816`*^9}, {3.7764093224678173`*^9, 
   3.776409572396212*^9}, {3.7764096224396653`*^9, 3.7764096638050823`*^9}, {
   3.7764097399995947`*^9, 3.7764097578963833`*^9}, {3.776409799262787*^9, 
   3.776409833496255*^9}, {3.7764099851317463`*^9, 3.77641000493445*^9}, {
   3.7764100382314663`*^9, 3.7764100440411367`*^9}, {3.7764101549478645`*^9, 
   3.7764101610623755`*^9}, {3.7764104033728766`*^9, 3.776410436862774*^9}, {
   3.7764104795384417`*^9, 3.776410549721387*^9}, 3.776428331968004*^9, {
   3.7764371243388977`*^9, 3.7764371303715067`*^9}, 3.776437164978725*^9, {
   3.77643721536313*^9, 3.7764372240422792`*^9}, {3.776439075019648*^9, 
   3.776439140184679*^9}, {3.776439208387807*^9, 3.776439214867412*^9}, {
   3.776439386368885*^9, 3.776439418350931*^9}, {3.776439468945496*^9, 
   3.776439476734679*^9}, {3.7764554244549875`*^9, 3.7764554634396553`*^9}, 
   3.7805092454427614`*^9, {3.7805092958725476`*^9, 3.780509303125936*^9}, {
   3.7805592832573133`*^9, 
   3.7805592908677907`*^9}},ExpressionUUID->"7e669c3b-775b-411b-b444-\
cc600741cad6"],

Cell[BoxData[
 TagBox[
  StyleBox[
   DynamicModuleBox[{$CellContext`a$$ = 1, $CellContext`b$$ = 
    0.365, $CellContext`\[Lambda]$$ = 0.1, $CellContext`\[CurlyRho]$$ = 360, 
    Typeset`show$$ = True, Typeset`bookmarkList$$ = {}, 
    Typeset`bookmarkMode$$ = "Menu", Typeset`animator$$, Typeset`animvar$$ = 
    1, Typeset`name$$ = "\"untitled\"", Typeset`specs$$ = {{{
       Hold[$CellContext`a$$], 1, "A"}, 1, 1.5}, {{
       Hold[$CellContext`b$$], 0.365, "\[Beta]"}, 0.1, 0.6}, {{
       Hold[$CellContext`\[Lambda]$$], 0.1, "\[CapitalLambda]"}, 0.1, 0.2}, {{
       Hold[$CellContext`\[CurlyRho]$$], 360, "\[Theta]0"}, 350, 370}}, 
    Typeset`size$$ = {360., {175., 180.}}, Typeset`update$$ = 0, 
    Typeset`initDone$$, Typeset`skipInitDone$$ = 
    False, $CellContext`a$18493$$ = 0, $CellContext`b$18494$$ = 
    0, $CellContext`\[Lambda]$18495$$ = 0, $CellContext`\[CurlyRho]$18496$$ = 
    0}, 
    DynamicBox[Manipulate`ManipulateBoxes[
     1, StandardForm, 
      "Variables" :> {$CellContext`a$$ = 1, $CellContext`b$$ = 
        0.365, $CellContext`\[Lambda]$$ = 0.1, $CellContext`\[CurlyRho]$$ = 
        360}, "ControllerVariables" :> {
        Hold[$CellContext`a$$, $CellContext`a$18493$$, 0], 
        Hold[$CellContext`b$$, $CellContext`b$18494$$, 0], 
        Hold[$CellContext`\[Lambda]$$, $CellContext`\[Lambda]$18495$$, 0], 
        Hold[$CellContext`\[CurlyRho]$$, $CellContext`\[CurlyRho]$18496$$, 
         0]}, "OtherVariables" :> {
       Typeset`show$$, Typeset`bookmarkList$$, Typeset`bookmarkMode$$, 
        Typeset`animator$$, Typeset`animvar$$, Typeset`name$$, 
        Typeset`specs$$, Typeset`size$$, Typeset`update$$, Typeset`initDone$$,
         Typeset`skipInitDone$$}, "Body" :> ContourPlot[
        ReplaceRepeated[
         Log[
          $CellContext`F[$CellContext`\[CapitalTheta], $CellContext`t]], \
{$CellContext`\[Theta]1 -> $CellContext`\[Theta]0 + \
$CellContext`A/$CellContext`\[CapitalLambda], $CellContext`\[Alpha] -> 
          1/
           13, $CellContext`\[Beta] -> $CellContext`b$$, $CellContext`A -> \
$CellContext`a$$, $CellContext`\[CapitalLambda] -> $CellContext`\[Lambda]$$, \
$CellContext`\[Theta]0 -> $CellContext`\[CurlyRho]$$}], {$CellContext`t, 0, 
         60}, {$CellContext`\[CapitalTheta], $CellContext`\[CurlyRho]$$, 
         400}], "Specifications" :> {{{$CellContext`a$$, 1, "A"}, 1, 1.5, 
         Appearance -> "Labeled"}, {{$CellContext`b$$, 0.365, "\[Beta]"}, 0.1,
          0.6, Appearance -> 
         "Labeled"}, {{$CellContext`\[Lambda]$$, 0.1, "\[CapitalLambda]"}, 
         0.1, 0.2, Appearance -> 
         "Labeled"}, {{$CellContext`\[CurlyRho]$$, 360, "\[Theta]0"}, 350, 
         370, Appearance -> "Labeled"}}, "Options" :> {}, 
      "DefaultOptions" :> {}],
     ImageSizeCache->{411., {264., 270.}},
     SingleEvaluation->True],
    Deinitialization:>None,
    DynamicModuleValues:>{},
    Initialization:>({$CellContext`F[
         Pattern[$CellContext`\[Theta], 
          Blank[]], 
         Pattern[$CellContext`t, 
          Blank[]]] := If[$CellContext`\[Theta] > $CellContext`\[Theta]1, 
         If[$CellContext`t > ($CellContext`\[Theta] - \
$CellContext`\[Theta]1)/$CellContext`A, 
          
          Exp[(-$CellContext`\[Alpha]) $CellContext`t - $CellContext`\
\[CapitalLambda] ($CellContext`t - ($CellContext`\[Theta] - $CellContext`\
\[Theta]1)/$CellContext`A) - (($CellContext`\[Beta] \
$CellContext`A)/$CellContext`\[CapitalLambda]) 
           Exp[(-$CellContext`\[CapitalLambda]) ($CellContext`t - \
($CellContext`\[Theta] - $CellContext`\[Theta]1)/$CellContext`A)]], 
          
          Exp[(-($CellContext`\[Alpha] - $CellContext`\[Beta] \
$CellContext`A)) $CellContext`t - $CellContext`\[Beta] ($CellContext`\[Theta] - \
$CellContext`\[Theta]0)]], 
         Exp[(-($CellContext`\[Alpha] + $CellContext`\[CapitalLambda])) \
$CellContext`t - ($CellContext`\[Beta] ($CellContext`\[Theta] - $CellContext`\
\[Theta]0)) Exp[(-$CellContext`\[CapitalLambda]) $CellContext`t]]]}; 
     Typeset`initDone$$ = True),
    SynchronousInitialization->True,
    UndoTrackedVariables:>{Typeset`show$$, Typeset`bookmarkMode$$},
    UnsavedVariables:>{Typeset`initDone$$},
    UntrackedVariables:>{Typeset`size$$}], "Manipulate",
   Deployed->True,
   StripOnInput->False],
  Manipulate`InterpretManipulate[1]]], "Output",
 CellChangeTimes->{3.776439429407762*^9, 3.7764394796800013`*^9, 
  3.7764539471129956`*^9, 3.7764554646419487`*^9, 3.7765080262144623`*^9, 
  3.7765812182862015`*^9, 3.7767066323078775`*^9, 3.7767071747953806`*^9, 
  3.776707401229658*^9, 3.780309165378862*^9, 3.7804223092171288`*^9, 
  3.7805093039383745`*^9, 3.780514428356386*^9, 3.7805149282125254`*^9, 
  3.780560760641322*^9},ExpressionUUID->"3cdf87bb-0244-4bb2-a9d9-\
6d02ad90c767"]
}, Open  ]]
}, Open  ]],

Cell[CellGroupData[{

Cell["Integral of the impact over time, normalization constant", "Chapter",
 CellChangeTimes->{{3.776411964643236*^9, 3.7764119905884347`*^9}, {
   3.7764284828879137`*^9, 3.776428500361361*^9}, {3.776428914916774*^9, 
   3.7764289253572474`*^9}, {3.776437680426756*^9, 3.776437693866437*^9}, 
   3.7764389439851418`*^9, 
   3.7764541379835625`*^9},ExpressionUUID->"65889081-75c0-4f05-ae9d-\
3ced8225cccc"],

Cell[CellGroupData[{

Cell["For \[Theta] < \[Theta]1", "Section",
 CellChangeTimes->{{3.776437700204116*^9, 3.776437701066161*^9}, {
   3.776437735580882*^9, 3.7764377556899548`*^9}, 
   3.7764541459180126`*^9},ExpressionUUID->"49a9ea81-f79e-4cbc-b27a-\
1cc93ffdabd7"],

Cell[CellGroupData[{

Cell[BoxData[
 RowBox[{"I0", " ", "=", " ", 
  RowBox[{"Integrate", "[", 
   RowBox[{
    RowBox[{"Exp", "[", 
     RowBox[{
      RowBox[{
       RowBox[{"-", 
        RowBox[{"(", 
         RowBox[{"\[Alpha]", "+", "\[CapitalLambda]"}], ")"}]}], "t"}], " ", 
      "-", " ", 
      RowBox[{"\[Beta]", 
       RowBox[{"(", 
        RowBox[{"\[Theta]", "-", "\[Theta]0"}], ")"}], " ", 
       RowBox[{"Exp", "[", 
        RowBox[{
         RowBox[{"-", "\[CapitalLambda]"}], " ", "t"}], "]"}]}]}], "]"}], ",", 
    RowBox[{"{", 
     RowBox[{"t", ",", "0", ",", "\[Infinity]"}], "}"}]}], "]"}]}]], "Input",
 CellChangeTimes->{{3.776437760166905*^9, 3.776437760441885*^9}, {
  3.776437794316434*^9, 3.776437828825622*^9}, {3.7804623755110474`*^9, 
  3.7804623784968495`*^9}},ExpressionUUID->"335a8e78-088d-4042-af62-\
4d1b067c96b7"],

Cell[BoxData[
 FractionBox[
  RowBox[{
   SuperscriptBox[
    RowBox[{"(", 
     RowBox[{"\[Beta]", " ", 
      RowBox[{"(", 
       RowBox[{"\[Theta]", "-", "\[Theta]0"}], ")"}]}], ")"}], 
    RowBox[{"-", 
     FractionBox[
      RowBox[{"\[Alpha]", "+", "\[CapitalLambda]"}], "\[CapitalLambda]"]}]], 
   " ", 
   RowBox[{"(", 
    RowBox[{
     RowBox[{"Gamma", "[", 
      FractionBox[
       RowBox[{"\[Alpha]", "+", "\[CapitalLambda]"}], "\[CapitalLambda]"], 
      "]"}], "-", 
     RowBox[{"Gamma", "[", 
      RowBox[{
       FractionBox[
        RowBox[{"\[Alpha]", "+", "\[CapitalLambda]"}], "\[CapitalLambda]"], 
       ",", 
       RowBox[{"\[Beta]", " ", 
        RowBox[{"(", 
         RowBox[{"\[Theta]", "-", "\[Theta]0"}], ")"}]}]}], "]"}]}], ")"}]}], 
  "\[CapitalLambda]"]], "Output",
 CellChangeTimes->{3.7764378524558268`*^9, 
  3.7804623922347794`*^9},ExpressionUUID->"0939f8c4-a8a7-48cd-b9d9-\
6e8d83ee3ba2"]
}, Open  ]]
}, Open  ]],

Cell[CellGroupData[{

Cell["For \[Theta] > \[Theta]1", "Section",
 CellChangeTimes->{{3.776437875226406*^9, 3.7764379444249496`*^9}, {
  3.776437975121181*^9, 3.776438099883463*^9}, {3.776454156469983*^9, 
  3.7764541666571245`*^9}},ExpressionUUID->"23015b51-d638-4ed6-bfbd-\
4f92d0bb9316"],

Cell[CellGroupData[{

Cell["Second part t > (\[Theta]-\[Theta]1)/A  (I2)", "Subsection",
 CellChangeTimes->{{3.776437875226406*^9, 3.7764379444249496`*^9}, {
  3.776437975121181*^9, 3.776438099883463*^9}, {3.776454156469983*^9, 
  3.77645417691824*^9}, {3.776454856220769*^9, 3.7764548619115057`*^9}, {
  3.7804614929054327`*^9, 
  3.7804615127043886`*^9}},ExpressionUUID->"c6ee4b5d-65c5-455b-8e78-\
dfd3ac080cfc"],

Cell["Integral after time shift to set the lower boundary to 0", "Text",
 CellChangeTimes->{{3.780318195604205*^9, 
  3.780318241246935*^9}},ExpressionUUID->"fdd5383f-495b-49e6-aa94-\
6c9e33ed4faa"],

Cell[CellGroupData[{

Cell[BoxData[
 RowBox[{"I2", " ", "=", " ", 
  RowBox[{
   RowBox[{"Exp", "[", 
    RowBox[{
     RowBox[{"-", "\[Alpha]"}], " ", 
     FractionBox[
      RowBox[{"\[Theta]", "-", "\[Theta]1"}], "A"]}], "]"}], 
   RowBox[{"Assuming", "[", 
    RowBox[{
     RowBox[{"{", 
      RowBox[{
       RowBox[{
        RowBox[{"\[Alpha]", "+", "\[CapitalLambda]"}], ">", "0"}], ",", 
       RowBox[{
        FractionBox[
         RowBox[{"\[Beta]", " ", "A"}], "\[CapitalLambda]"], ">", "0"}], ",", 
       RowBox[{"\[CapitalLambda]", ">", "0"}]}], "}"}], ",", 
     RowBox[{"Integrate", "[", 
      RowBox[{
       RowBox[{"Exp", "[", 
        RowBox[{
         RowBox[{
          RowBox[{"-", 
           RowBox[{"(", 
            RowBox[{"\[Alpha]", "+", "\[CapitalLambda]"}], ")"}]}], " ", 
          "t"}], " ", "-", 
         RowBox[{
          FractionBox[
           RowBox[{"\[Beta]", " ", "A"}], "\[CapitalLambda]"], 
          RowBox[{"Exp", "[", 
           RowBox[{
            RowBox[{"-", "\[CapitalLambda]"}], " ", "t"}], "]"}]}]}], "]"}], 
       ",", 
       RowBox[{"{", 
        RowBox[{"t", ",", "0", ",", "\[Infinity]"}], "}"}]}], "]"}]}], 
    "]"}]}]}]], "Input",
 CellChangeTimes->{{3.776437875226406*^9, 3.7764379444249496`*^9}, {
  3.776437975121181*^9, 3.7764381614977713`*^9}, {3.776438242801354*^9, 
  3.7764382662390337`*^9}, {3.776438310448892*^9, 3.776438319478526*^9}, {
  3.780461477262925*^9, 3.7804614809959087`*^9}, {3.780461984728771*^9, 
  3.780461987464411*^9}, {3.780463015019453*^9, 
  3.7804630159881763`*^9}},ExpressionUUID->"45788ab6-b913-4357-bb5e-\
93b46866de11"],

Cell[BoxData[
 RowBox[{
  SuperscriptBox["\[ExponentialE]", 
   RowBox[{"-", 
    FractionBox[
     RowBox[{"\[Alpha]", " ", 
      RowBox[{"(", 
       RowBox[{"\[Theta]", "-", "\[Theta]1"}], ")"}]}], "A"]}]], " ", 
  SuperscriptBox[
   RowBox[{"(", 
    RowBox[{"A", " ", "\[Beta]"}], ")"}], 
   RowBox[{"-", 
    FractionBox[
     RowBox[{"\[Alpha]", "+", "\[CapitalLambda]"}], "\[CapitalLambda]"]}]], 
  " ", 
  SuperscriptBox["\[CapitalLambda]", 
   RowBox[{"\[Alpha]", "/", "\[CapitalLambda]"}]], " ", 
  RowBox[{"(", 
   RowBox[{
    RowBox[{"Gamma", "[", 
     FractionBox[
      RowBox[{"\[Alpha]", "+", "\[CapitalLambda]"}], "\[CapitalLambda]"], 
     "]"}], "-", 
    RowBox[{"Gamma", "[", 
     RowBox[{
      FractionBox[
       RowBox[{"\[Alpha]", "+", "\[CapitalLambda]"}], "\[CapitalLambda]"], 
      ",", 
      FractionBox[
       RowBox[{"A", " ", "\[Beta]"}], "\[CapitalLambda]"]}], "]"}]}], 
   ")"}]}]], "Output",
 CellChangeTimes->{{3.776438102229121*^9, 3.7764381104128733`*^9}, 
   3.7764381539172773`*^9, 3.776438214141924*^9, 3.776438292598301*^9, 
   3.776438324500046*^9, 3.7804620113934593`*^9, 3.780463031305398*^9, 
   3.780514504648638*^9},ExpressionUUID->"0a2ad285-f54e-4da7-afb1-\
c74f89f642d4"]
}, Open  ]],

Cell[TextData[{
 "Exponential decay as a function of \[Theta] as ",
 Cell[BoxData[
  SuperscriptBox["\[ExponentialE]", 
   RowBox[{"-", 
    FractionBox[
     RowBox[{"\[Alpha]", " ", 
      RowBox[{"(", 
       RowBox[{"\[Theta]", "-", "\[Theta]1"}], ")"}]}], "A"]}]]],
  CellChangeTimes->{{3.776438102229121*^9, 3.7764381104128733`*^9}, 
    3.7764381539172773`*^9, 3.776438214141924*^9, 3.776438292598301*^9, 
    3.776438324500046*^9, 3.7804620113934593`*^9, 3.780463031305398*^9},
  ExpressionUUID->"36d60791-9890-4a59-807f-7d6e05bfe668"]
}], "Text",
 CellChangeTimes->{{3.780487038755517*^9, 
  3.7804871015174637`*^9}},ExpressionUUID->"50496439-cb9f-48c2-8f54-\
325085d4afa0"]
}, Open  ]],

Cell[CellGroupData[{

Cell["First part t < (\[Theta]-\[Theta]1)/A  (I1)", "Subsection",
 CellChangeTimes->{{3.7764381974155493`*^9, 3.7764382285348377`*^9}, {
   3.7764383491095667`*^9, 3.77643836104848*^9}, 3.7764541851914964`*^9, {
   3.776454737663742*^9, 3.7764547475700808`*^9}, {3.7764548317607927`*^9, 
   3.77645484913785*^9}, {3.7804615009855795`*^9, 
   3.78046151889258*^9}},ExpressionUUID->"4f709d11-8991-4045-8acf-\
83ddc3a763f6"],

Cell[CellGroupData[{

Cell[BoxData[
 RowBox[{"I1", " ", "=", " ", 
  RowBox[{"FullSimplify", "[", 
   RowBox[{"Integrate", "[", 
    RowBox[{
     RowBox[{"Exp", "[", 
      RowBox[{
       RowBox[{
        RowBox[{"-", 
         RowBox[{"(", 
          RowBox[{"\[Alpha]", " ", "-", " ", 
           RowBox[{"\[Beta]", " ", "A"}]}], ")"}]}], "t"}], " ", "-", 
       RowBox[{"\[Beta]", 
        RowBox[{"(", 
         RowBox[{"\[Theta]", "-", "\[Theta]0"}], ")"}]}]}], "]"}], ",", 
     RowBox[{"{", 
      RowBox[{"t", ",", "0", ",", 
       FractionBox[
        RowBox[{"\[Theta]", "-", "\[Theta]1"}], "A"]}], "}"}]}], "]"}], " ", 
   "]"}]}]], "Input",
 CellChangeTimes->{{3.7764381974155493`*^9, 3.7764382285348377`*^9}, {
  3.7764383491095667`*^9, 3.7764383643333673`*^9}, {3.776438402760672*^9, 
  3.776438464172987*^9}, {3.776454257112261*^9, 3.77645427743363*^9}, {
  3.7764543106015935`*^9, 3.776454321577301*^9}, {3.7804618196394854`*^9, 
  3.780461852522418*^9}, {3.7804619932296295`*^9, 3.780461995932499*^9}, {
  3.7804630194583044`*^9, 
  3.7804630231641083`*^9}},ExpressionUUID->"fedcd1fb-8b64-467e-8ec0-\
001ed6ab0fa0"],

Cell[BoxData[
 FractionBox[
  RowBox[{
   SuperscriptBox["\[ExponentialE]", 
    RowBox[{"\[Beta]", " ", 
     RowBox[{"(", 
      RowBox[{
       RowBox[{"-", "\[Theta]"}], "+", "\[Theta]0"}], ")"}]}]], " ", 
   RowBox[{"(", 
    RowBox[{
     RowBox[{"-", "1"}], "+", 
     SuperscriptBox["\[ExponentialE]", 
      FractionBox[
       RowBox[{
        RowBox[{"(", 
         RowBox[{
          RowBox[{"-", "\[Alpha]"}], "+", 
          RowBox[{"A", " ", "\[Beta]"}]}], ")"}], " ", 
        RowBox[{"(", 
         RowBox[{"\[Theta]", "-", "\[Theta]1"}], ")"}]}], "A"]]}], ")"}]}], 
  RowBox[{
   RowBox[{"-", "\[Alpha]"}], "+", 
   RowBox[{"A", " ", "\[Beta]"}]}]]], "Output",
 CellChangeTimes->{
  3.780514499742738*^9},ExpressionUUID->"7d22e2b1-8524-45e3-8424-\
0ea61d1181a5"]
}, Open  ]],

Cell[TextData[{
 "Notice that I1 = 0 for \[Theta] = \[Theta]1 and is always positive for \
\[Theta] > \[Theta]1\nIf \[Beta] > ",
 Cell[BoxData[
  FormBox[
   FractionBox["\[Alpha]", "A"], TraditionalForm]],ExpressionUUID->
  "ec441c8b-9d03-42de-805b-f6930c2ae731"],
 ", decays as  ",
 Cell[BoxData[
  SuperscriptBox["\[ExponentialE]", 
   RowBox[{"-", 
    FractionBox[
     RowBox[{"\[Alpha]", " ", 
      RowBox[{"(", 
       RowBox[{"\[Theta]", "-", "\[Theta]1"}], ")"}]}], "A"]}]]],
  CellChangeTimes->{{3.776438102229121*^9, 3.7764381104128733`*^9}, 
    3.7764381539172773`*^9, 3.776438214141924*^9, 3.776438292598301*^9, 
    3.776438324500046*^9, 3.7804620113934593`*^9, 3.780463031305398*^9},
  ExpressionUUID->"0b67da3a-7091-4181-bef9-f31dbcd0671d"],
 "at large \[Theta]. The ratio I1/I2 therefore tends to a constant"
}], "Text",
 CellChangeTimes->{{3.7804808878792725`*^9, 3.780480890935537*^9}, {
   3.780481044669242*^9, 3.780481169501337*^9}, 3.7804836861657295`*^9, {
   3.7804871990334673`*^9, 3.7804872010623136`*^9}, {3.7804875667103033`*^9, 
   3.780487700213571*^9}, {3.7804877821624975`*^9, 
   3.780487892328257*^9}},ExpressionUUID->"784746ff-374c-44af-94b6-\
8b0ad1f03400"]
}, Open  ]]
}, Open  ]],

Cell[CellGroupData[{

Cell["Generic integral and rewritting of I0 and I2", "Section",
 CellChangeTimes->{{3.776454193477767*^9, 3.7764542061774645`*^9}, {
  3.7804624385035424`*^9, 3.7804624444122972`*^9}, {3.7804626516731896`*^9, 
  3.7804626569853125`*^9}, {3.78046341199571*^9, 
  3.780463414015959*^9}},ExpressionUUID->"869de240-e246-4f88-8a84-\
8b69843583c3"],

Cell[CellGroupData[{

Cell[BoxData[
 RowBox[{
  RowBox[{"G", "[", 
   RowBox[{"a_", ",", "\[Xi]_"}], "]"}], "=", 
  RowBox[{"Evaluate", "[", 
   RowBox[{"Assuming", "[", 
    RowBox[{
     RowBox[{"a", ">", "0"}], ",", 
     RowBox[{"Integrate", "[", 
      RowBox[{
       RowBox[{"Exp", "[", 
        RowBox[{
         RowBox[{
          RowBox[{"-", " ", "a"}], " ", "t"}], " ", "-", 
         RowBox[{"\[Xi]", " ", 
          RowBox[{"Exp", "[", 
           RowBox[{"-", "t"}], "]"}]}]}], "]"}], ",", 
       RowBox[{"{", 
        RowBox[{"t", ",", "0", ",", "\[Infinity]"}], "}"}]}], "]"}]}], "]"}], 
   "]"}]}]], "Input",
 CellChangeTimes->{{3.7764539371307197`*^9, 3.7764539405047884`*^9}, {
  3.7803364231617193`*^9, 3.780336440407691*^9}, {3.7803364705960555`*^9, 
  3.7803365637624063`*^9}},ExpressionUUID->"6a8404d0-4df9-4fda-8b62-\
338f5ac43052"],

Cell[BoxData[
 RowBox[{
  SuperscriptBox["\[Xi]", 
   RowBox[{"-", "a"}]], " ", 
  RowBox[{"(", 
   RowBox[{
    RowBox[{"Gamma", "[", "a", "]"}], "-", 
    RowBox[{"Gamma", "[", 
     RowBox[{"a", ",", "\[Xi]"}], "]"}]}], ")"}]}]], "Output",
 CellChangeTimes->{3.7803365666949615`*^9, 3.780464805038591*^9, 
  3.7805144322623606`*^9, 3.7805149328371954`*^9, 
  3.780559452243682*^9},ExpressionUUID->"5fa078d3-7cfd-4d73-b13e-\
38dae4e7acc2"]
}, Open  ]],

Cell[CellGroupData[{

Cell["\<\
I1 and G are repeated to avoid calculations in the initialization\
\>", "Subsection",
 CellChangeTimes->{{3.780559528312146*^9, 
  3.7805595669524727`*^9}},ExpressionUUID->"8ddc7c60-240a-4a87-8251-\
9a209e283a4a"],

Cell[BoxData[{
 RowBox[{
  RowBox[{"I0", " ", "=", " ", 
   RowBox[{
    FractionBox["1", "\[CapitalLambda]"], 
    RowBox[{"G", "[", 
     RowBox[{
      FractionBox[
       RowBox[{"\[Alpha]", "+", "\[CapitalLambda]"}], "\[CapitalLambda]"], 
      ",", 
      RowBox[{"\[Beta]", " ", 
       RowBox[{"(", 
        RowBox[{"\[Theta]", "-", "\[Theta]0"}], ")"}]}]}], "]"}]}]}], 
  ";"}], "\[IndentingNewLine]", 
 RowBox[{
  RowBox[{"I2", " ", "=", " ", 
   RowBox[{
    SuperscriptBox["\[ExponentialE]", 
     RowBox[{"-", 
      FractionBox[
       RowBox[{"\[Alpha]", " ", 
        RowBox[{"(", 
         RowBox[{"\[Theta]", "-", "\[Theta]1"}], ")"}]}], "A"]}]], 
    FractionBox["1", "\[CapitalLambda]"], 
    RowBox[{"G", "[", 
     RowBox[{
      FractionBox[
       RowBox[{"\[Alpha]", "+", "\[CapitalLambda]"}], "\[CapitalLambda]"], 
      ",", 
      FractionBox[
       RowBox[{"A", " ", "\[Beta]"}], "\[CapitalLambda]"]}], "]"}]}]}], 
  ";"}], "\[IndentingNewLine]", 
 RowBox[{
  RowBox[{"I1", " ", "=", " ", 
   FractionBox[
    RowBox[{
     SuperscriptBox["\[ExponentialE]", 
      RowBox[{"\[Beta]", " ", 
       RowBox[{"(", 
        RowBox[{
         RowBox[{"-", "\[Theta]"}], "+", "\[Theta]0"}], ")"}]}]], " ", 
     RowBox[{"(", 
      RowBox[{
       RowBox[{"-", "1"}], "+", 
       SuperscriptBox["\[ExponentialE]", 
        FractionBox[
         RowBox[{
          RowBox[{"(", 
           RowBox[{
            RowBox[{"-", "\[Alpha]"}], "+", 
            RowBox[{"A", " ", "\[Beta]"}]}], ")"}], " ", 
          RowBox[{"(", 
           RowBox[{"\[Theta]", "-", "\[Theta]1"}], ")"}]}], "A"]]}], ")"}]}], 
    RowBox[{
     RowBox[{"-", "\[Alpha]"}], "+", 
     RowBox[{"A", " ", "\[Beta]"}]}]]}], ";"}], "\[IndentingNewLine]", 
 RowBox[{
  RowBox[{
   RowBox[{"G", "[", 
    RowBox[{"a_", ",", "\[Xi]_"}], "]"}], ":=", " ", 
   RowBox[{
    SuperscriptBox["\[Xi]", 
     RowBox[{"-", "a"}]], " ", 
    RowBox[{"(", 
     RowBox[{
      RowBox[{"Gamma", "[", "a", "]"}], "-", 
      RowBox[{"Gamma", "[", 
       RowBox[{"a", ",", "\[Xi]"}], "]"}]}], ")"}]}]}], ";"}]}], "Input",
 InitializationCell->True,
 CellChangeTimes->{{3.780462459960874*^9, 3.7804625746411257`*^9}, {
   3.780462604871255*^9, 3.780462624264604*^9}, {3.7804626796923428`*^9, 
   3.780462680754768*^9}, 3.7804631129996295`*^9, {3.780463417609451*^9, 
   3.7804634183905926`*^9}, {3.780559332362756*^9, 3.7805593447296944`*^9}, 
   3.780559455868498*^9, {3.7805594972648573`*^9, 3.7805595100007334`*^9}, {
   3.7805597247054157`*^9, 
   3.780559725924145*^9}},ExpressionUUID->"c7696fcc-cdb5-4ba2-b455-\
2b9760e13587"],

Cell[BoxData[{
 RowBox[{
  RowBox[{"Total", " ", "for", " ", "\[Theta]"}], ">", 
  RowBox[{"\[Theta]1", " ", 
   RowBox[{"(", 
    RowBox[{"I1", " ", "+", "I2"}], ")"}]}]}], "\[IndentingNewLine]", 
 RowBox[{" ", 
  RowBox[{"n2", "=", 
   RowBox[{
    FractionBox[
     RowBox[{
      SuperscriptBox["\[ExponentialE]", 
       RowBox[{"\[Beta]", " ", 
        RowBox[{"(", 
         RowBox[{
          RowBox[{"-", "\[Theta]"}], "+", "\[Theta]0"}], ")"}]}]], " ", 
      RowBox[{"(", 
       RowBox[{
        RowBox[{"-", "1"}], "+", 
        SuperscriptBox["\[ExponentialE]", 
         FractionBox[
          RowBox[{
           RowBox[{"(", 
            RowBox[{
             RowBox[{"-", "\[Alpha]"}], "+", 
             RowBox[{"A", " ", "\[Beta]"}]}], ")"}], " ", 
           RowBox[{"(", 
            RowBox[{"\[Theta]", "-", "\[Theta]1"}], ")"}]}], "A"]]}], ")"}]}], 
     RowBox[{
      RowBox[{"-", "\[Alpha]"}], "+", 
      RowBox[{"A", " ", "\[Beta]"}]}]], "+", 
    RowBox[{
     SuperscriptBox["\[ExponentialE]", 
      RowBox[{"-", 
       FractionBox[
        RowBox[{"\[Alpha]", " ", 
         RowBox[{"(", 
          RowBox[{"\[Theta]", "-", "\[Theta]1"}], ")"}]}], "A"]}]], " ", 
     FractionBox["1", "\[CapitalLambda]"], 
     RowBox[{"G", "[", 
      RowBox[{
       FractionBox[
        RowBox[{"\[Alpha]", "+", "\[CapitalLambda]"}], "\[CapitalLambda]"], 
       ",", 
       FractionBox[
        RowBox[{"A", " ", "\[Beta]"}], "\[CapitalLambda]"]}], 
      "]"}]}]}]}]}]}], "Text",
 CellChangeTimes->{{3.776454872195611*^9, 3.7764549272320557`*^9}, {
   3.7764550311055107`*^9, 3.776455115208311*^9}, {3.776706787179356*^9, 
   3.7767067892091913`*^9}, 3.780336379276536*^9, {3.780461441922532*^9, 
   3.7804614515312786`*^9}, {3.7804615213923535`*^9, 3.780461533403143*^9}, {
   3.7804618643444386`*^9, 3.7804618720171385`*^9}, {3.7804631047773857`*^9, 
   3.7804631307221136`*^9}},ExpressionUUID->"86ce030f-ef92-4766-84a1-\
fc817d163ccf"]
}, Open  ]]
}, Open  ]]
}, Open  ]],

Cell[CellGroupData[{

Cell["Normalized impact (dividing by the integral in t)", "Chapter",
 CellChangeTimes->{{3.7764540436326585`*^9, 3.7764540946654167`*^9}, 
   3.776454227316355*^9},ExpressionUUID->"7fad5e09-8f8e-474b-9c45-\
f7396711b77f"],

Cell[CellGroupData[{

Cell["Normalized impact formula", "Section",
 CellChangeTimes->{{3.776516593843047*^9, 
  3.7765166091263027`*^9}},ExpressionUUID->"5efd3c13-e37e-4e02-9547-\
65f44ca71744"],

Cell[BoxData[
 RowBox[{
  RowBox[{"Fn", "[", 
   RowBox[{"\[Theta]_", ",", "t_"}], "]"}], " ", ":=", " ", 
  RowBox[{"If", "[", 
   RowBox[{
    RowBox[{"\[Theta]", ">", "\[Theta]1"}], ",", 
    RowBox[{
     RowBox[{"If", "[", 
      RowBox[{
       RowBox[{"t", ">", 
        FractionBox[
         RowBox[{"\[Theta]", "-", "\[Theta]1"}], "A"]}], ",", 
       RowBox[{"Exp", "[", 
        RowBox[{
         RowBox[{
          RowBox[{"-", "\[Alpha]"}], " ", "t"}], "-", 
         RowBox[{"\[CapitalLambda]", 
          RowBox[{"(", 
           RowBox[{"t", "-", 
            FractionBox[
             RowBox[{"\[Theta]", "-", "\[Theta]1"}], "A"]}], ")"}]}], " ", 
         "-", 
         RowBox[{
          FractionBox[
           RowBox[{"\[Beta]", " ", "A"}], "\[CapitalLambda]"], 
          RowBox[{"Exp", "[", 
           RowBox[{
            RowBox[{"-", "\[CapitalLambda]"}], 
            RowBox[{"(", 
             RowBox[{"t", "-", 
              FractionBox[
               RowBox[{"\[Theta]", "-", "\[Theta]1"}], "A"]}], ")"}]}], 
           "]"}]}]}], "]"}], ",", 
       RowBox[{"Exp", "[", 
        RowBox[{
         RowBox[{
          RowBox[{"-", 
           RowBox[{"(", 
            RowBox[{"\[Alpha]", " ", "-", " ", 
             RowBox[{"\[Beta]", " ", "A"}]}], ")"}]}], "t"}], " ", "-", 
         RowBox[{"\[Beta]", 
          RowBox[{"(", 
           RowBox[{"\[Theta]", "-", "\[Theta]0"}], ")"}]}]}], "]"}]}], "]"}], 
     "/", " ", 
     RowBox[{"(", 
      RowBox[{
       FractionBox[
        RowBox[{
         SuperscriptBox["\[ExponentialE]", 
          RowBox[{"\[Beta]", " ", 
           RowBox[{"(", 
            RowBox[{
             RowBox[{"-", "\[Theta]"}], "+", "\[Theta]0"}], ")"}]}]], " ", 
         RowBox[{"(", 
          RowBox[{
           RowBox[{"-", "1"}], "+", 
           SuperscriptBox["\[ExponentialE]", 
            FractionBox[
             RowBox[{
              RowBox[{"(", 
               RowBox[{
                RowBox[{"-", "\[Alpha]"}], "+", 
                RowBox[{"A", " ", "\[Beta]"}]}], ")"}], " ", 
              RowBox[{"(", 
               RowBox[{"\[Theta]", "-", "\[Theta]1"}], ")"}]}], "A"]]}], 
          ")"}]}], 
        RowBox[{
         RowBox[{"-", "\[Alpha]"}], "+", 
         RowBox[{"A", " ", "\[Beta]"}]}]], "+", 
       RowBox[{
        SuperscriptBox["\[ExponentialE]", 
         RowBox[{"-", 
          FractionBox[
           RowBox[{"\[Alpha]", " ", 
            RowBox[{"(", 
             RowBox[{"\[Theta]", "-", "\[Theta]1"}], ")"}]}], "A"]}]], 
        FractionBox["1", "\[CapitalLambda]"], 
        RowBox[{"G", "[", 
         RowBox[{
          FractionBox[
           RowBox[{"\[Alpha]", "+", "\[CapitalLambda]"}], "\[CapitalLambda]"],
           ",", 
          FractionBox[
           RowBox[{"A", " ", "\[Beta]"}], "\[CapitalLambda]"]}], "]"}]}]}], 
      ")"}]}], ",", 
    FractionBox[
     RowBox[{"\[CapitalLambda]", " ", 
      RowBox[{"Exp", "[", 
       RowBox[{
        RowBox[{
         RowBox[{"-", 
          RowBox[{"(", 
           RowBox[{"\[Alpha]", "+", "\[CapitalLambda]"}], ")"}]}], "t"}], " ",
         "-", " ", 
        RowBox[{"\[Beta]", 
         RowBox[{"(", 
          RowBox[{"\[Theta]", "-", "\[Theta]0"}], ")"}], " ", 
         RowBox[{"Exp", "[", 
          RowBox[{
           RowBox[{"-", "\[CapitalLambda]"}], " ", "t"}], "]"}]}]}], "]"}]}], 
     RowBox[{"G", "[", 
      RowBox[{
       FractionBox[
        RowBox[{"\[Alpha]", "+", "\[CapitalLambda]"}], "\[CapitalLambda]"], 
       ",", 
       RowBox[{"\[Beta]", " ", 
        RowBox[{"(", 
         RowBox[{"\[Theta]", "-", "\[Theta]0"}], ")"}]}]}], "]"}]]}], 
   "]"}]}]], "Input",
 InitializationCell->True,
 CellChangeTimes->{{3.7764540039733953`*^9, 3.776454012188701*^9}, {
   3.776454112589122*^9, 3.7764541198639684`*^9}, {3.7764543534660172`*^9, 
   3.7764543633683395`*^9}, {3.7764544007599015`*^9, 
   3.7764544070852757`*^9}, {3.7764545317388096`*^9, 
   3.7764545380152082`*^9}, {3.776454581603217*^9, 3.7764546410451374`*^9}, 
   3.776454673200713*^9, {3.7764552023983126`*^9, 3.776455224228788*^9}, {
   3.7764592379356847`*^9, 3.776459344737432*^9}, {3.7765181577527103`*^9, 
   3.7765181779201756`*^9}, {3.7803367451774807`*^9, 3.780336753021762*^9}, {
   3.7803374850599256`*^9, 3.780337496229148*^9}, 
   3.7804620435787554`*^9},ExpressionUUID->"62a3ecfb-ddaa-40c4-bfaf-\
e3a256be4af9"]
}, Open  ]],

Cell[CellGroupData[{

Cell["Age / altitude map", "Section",
 CellChangeTimes->{{3.776516632774742*^9, 
  3.776516673772273*^9}},ExpressionUUID->"2c660588-a942-43e3-a1da-\
1b64c15827aa"],

Cell[CellGroupData[{

Cell[BoxData[
 RowBox[{"Manipulate", "[", 
  RowBox[{
   RowBox[{"ContourPlot", "[", 
    RowBox[{
     RowBox[{
      RowBox[{"Log", "[", 
       RowBox[{"Fn", "[", 
        RowBox[{"\[CapitalTheta]", ",", "t"}], "]"}], "]"}], "//.", 
      "\[VeryThinSpace]", 
      RowBox[{"{", "\[VeryThinSpace]", 
       RowBox[{
        RowBox[{"\[Theta]1", "\[Rule]", 
         RowBox[{"\[Theta]0", "+", " ", 
          FractionBox["A", "\[CapitalLambda]"]}]}], ",", 
        RowBox[{"\[Alpha]", "\[Rule]", 
         FractionBox["1", "13"]}], ",", 
        RowBox[{"\[Beta]", "\[Rule]", "b"}], ",", 
        RowBox[{"A", "\[Rule]", "a"}], ",", 
        RowBox[{"\[CapitalLambda]", "\[Rule]", "\[Lambda]"}], ",", 
        RowBox[{"\[Theta]0", "\[Rule]", "\[CurlyRho]"}]}], "}"}]}], ",", 
     RowBox[{"{", 
      RowBox[{"t", ",", "0", ",", "60"}], "}"}], ",", 
     RowBox[{"{", 
      RowBox[{"\[CapitalTheta]", ",", "\[CurlyRho]", ",", "440"}], "}"}]}], 
    "]"}], ",", 
   RowBox[{"{", 
    RowBox[{
     RowBox[{"{", 
      RowBox[{"a", ",", "1", ",", "\"\<A\>\""}], "}"}], ",", "1", ",", "1.5", 
     ",", 
     RowBox[{"Appearance", "\[Rule]", "\"\<Labeled\>\""}]}], "}"}], ",", 
   RowBox[{"{", 
    RowBox[{
     RowBox[{"{", 
      RowBox[{"b", ",", "0.365", ",", "\"\<\[Beta]\>\""}], "}"}], ",", "0.1", 
     ",", "0.6", ",", 
     RowBox[{"Appearance", "\[Rule]", "\"\<Labeled\>\""}]}], "}"}], ",", 
   RowBox[{"{", 
    RowBox[{
     RowBox[{"{", 
      RowBox[{"\[Lambda]", ",", "0.1", ",", "\"\<\[CapitalLambda]\>\""}], 
      "}"}], ",", "0.1", ",", "0.2", ",", 
     RowBox[{"Appearance", "\[Rule]", "\"\<Labeled\>\""}]}], "}"}], ",", 
   RowBox[{"{", 
    RowBox[{
     RowBox[{"{", 
      RowBox[{"\[CurlyRho]", ",", "360", ",", "\"\<\[Theta]0\>\""}], "}"}], 
     ",", "350", ",", "370", ",", 
     RowBox[{"Appearance", "\[Rule]", "\"\<Labeled\>\""}]}], "}"}]}], 
  "]"}]], "Input",
 InitializationCell->True,
 CellChangeTimes->{{3.7764089115682364`*^9, 3.7764090226098704`*^9}, {
   3.776409067331358*^9, 3.7764090779762816`*^9}, {3.7764093224678173`*^9, 
   3.776409572396212*^9}, {3.7764096224396653`*^9, 3.7764096638050823`*^9}, {
   3.7764097399995947`*^9, 3.7764097578963833`*^9}, {3.776409799262787*^9, 
   3.776409833496255*^9}, {3.7764099851317463`*^9, 3.77641000493445*^9}, {
   3.7764100382314663`*^9, 3.7764100440411367`*^9}, {3.7764101549478645`*^9, 
   3.7764101610623755`*^9}, {3.7764104033728766`*^9, 3.776410436862774*^9}, {
   3.7764104795384417`*^9, 3.776410549721387*^9}, 3.776428331968004*^9, {
   3.7764371243388977`*^9, 3.7764371303715067`*^9}, 3.776437164978725*^9, {
   3.77643721536313*^9, 3.7764372240422792`*^9}, {3.776439075019648*^9, 
   3.776439140184679*^9}, {3.776439208387807*^9, 3.776439214867412*^9}, {
   3.776439386368885*^9, 3.776439418350931*^9}, {3.776439468945496*^9, 
   3.776439476734679*^9}, {3.776455247384509*^9, 3.7764552752285447`*^9}, {
   3.7764593669127183`*^9, 3.776459380258067*^9}, {3.776459430257401*^9, 
   3.7764594310699353`*^9}, {3.776459507482125*^9, 3.776459513058928*^9}, 
   3.78033745642294*^9, {3.780337581852612*^9, 
   3.780337590710145*^9}},ExpressionUUID->"554bcbeb-5eca-435c-be63-\
b02d5784b0f6"],

Cell[BoxData[
 TagBox[
  StyleBox[
   DynamicModuleBox[{$CellContext`a$$ = 1, $CellContext`b$$ = 
    0.365, $CellContext`\[Lambda]$$ = 0.1, $CellContext`\[CurlyRho]$$ = 360, 
    Typeset`show$$ = True, Typeset`bookmarkList$$ = {}, 
    Typeset`bookmarkMode$$ = "Menu", Typeset`animator$$, Typeset`animvar$$ = 
    1, Typeset`name$$ = "\"untitled\"", Typeset`specs$$ = {{{
       Hold[$CellContext`a$$], 1, "A"}, 1, 1.5}, {{
       Hold[$CellContext`b$$], 0.365, "\[Beta]"}, 0.1, 0.6}, {{
       Hold[$CellContext`\[Lambda]$$], 0.1, "\[CapitalLambda]"}, 0.1, 0.2}, {{
       Hold[$CellContext`\[CurlyRho]$$], 360, "\[Theta]0"}, 350, 370}}, 
    Typeset`size$$ = {360., {175., 180.}}, Typeset`update$$ = 0, 
    Typeset`initDone$$, Typeset`skipInitDone$$ = 
    True, $CellContext`a$18562$$ = 0, $CellContext`b$18563$$ = 
    0, $CellContext`\[Lambda]$18564$$ = 0, $CellContext`\[CurlyRho]$18565$$ = 
    0}, 
    DynamicBox[Manipulate`ManipulateBoxes[
     1, StandardForm, 
      "Variables" :> {$CellContext`a$$ = 1, $CellContext`b$$ = 
        0.365, $CellContext`\[Lambda]$$ = 0.1, $CellContext`\[CurlyRho]$$ = 
        360}, "ControllerVariables" :> {
        Hold[$CellContext`a$$, $CellContext`a$18562$$, 0], 
        Hold[$CellContext`b$$, $CellContext`b$18563$$, 0], 
        Hold[$CellContext`\[Lambda]$$, $CellContext`\[Lambda]$18564$$, 0], 
        Hold[$CellContext`\[CurlyRho]$$, $CellContext`\[CurlyRho]$18565$$, 
         0]}, "OtherVariables" :> {
       Typeset`show$$, Typeset`bookmarkList$$, Typeset`bookmarkMode$$, 
        Typeset`animator$$, Typeset`animvar$$, Typeset`name$$, 
        Typeset`specs$$, Typeset`size$$, Typeset`update$$, Typeset`initDone$$,
         Typeset`skipInitDone$$}, "Body" :> ContourPlot[
        ReplaceRepeated[
         Log[
          $CellContext`Fn[$CellContext`\[CapitalTheta], $CellContext`t]], \
{$CellContext`\[Theta]1 -> $CellContext`\[Theta]0 + \
$CellContext`A/$CellContext`\[CapitalLambda], $CellContext`\[Alpha] -> 
          1/
           13, $CellContext`\[Beta] -> $CellContext`b$$, $CellContext`A -> \
$CellContext`a$$, $CellContext`\[CapitalLambda] -> $CellContext`\[Lambda]$$, \
$CellContext`\[Theta]0 -> $CellContext`\[CurlyRho]$$}], {$CellContext`t, 0, 
         60}, {$CellContext`\[CapitalTheta], $CellContext`\[CurlyRho]$$, 
         440}], "Specifications" :> {{{$CellContext`a$$, 1, "A"}, 1, 1.5, 
         Appearance -> "Labeled"}, {{$CellContext`b$$, 0.365, "\[Beta]"}, 0.1,
          0.6, Appearance -> 
         "Labeled"}, {{$CellContext`\[Lambda]$$, 0.1, "\[CapitalLambda]"}, 
         0.1, 0.2, Appearance -> 
         "Labeled"}, {{$CellContext`\[CurlyRho]$$, 360, "\[Theta]0"}, 350, 
         370, Appearance -> "Labeled"}}, "Options" :> {}, 
      "DefaultOptions" :> {}],
     ImageSizeCache->{411., {264., 270.}},
     SingleEvaluation->True],
    Deinitialization:>None,
    DynamicModuleValues:>{},
    SynchronousInitialization->True,
    UndoTrackedVariables:>{Typeset`show$$, Typeset`bookmarkMode$$},
    UnsavedVariables:>{Typeset`initDone$$},
    UntrackedVariables:>{Typeset`size$$}], "Manipulate",
   Deployed->True,
   StripOnInput->False],
  Manipulate`InterpretManipulate[1]]], "Output",
 CellChangeTimes->{
  3.776455276592764*^9, 3.7764592663843575`*^9, 3.776459382204971*^9, 
   3.7764594378770313`*^9, 3.776459514953843*^9, 3.776508027554694*^9, 
   3.7765083406783924`*^9, 3.776581219579462*^9, {3.7765814249121137`*^9, 
   3.7765814354211063`*^9}, 3.7767066343607044`*^9, 3.7767071749303007`*^9, 
   3.7767074013605633`*^9, 3.780309166660038*^9, {3.7803369413081846`*^9, 
   3.780336953928795*^9}, 3.7803369898005047`*^9, 3.7803374594305*^9, 
   3.78033759342864*^9, 3.780422310467039*^9, 3.780487929735774*^9, 
   3.7805144324810925`*^9, 3.780514934392868*^9, 3.78055939743976*^9, 
   3.7805607610631666`*^9},ExpressionUUID->"c1656431-dca3-493f-a989-\
bd33d7f5197a"]
}, Open  ]]
}, Open  ]],

Cell[CellGroupData[{

Cell["Section in t (age) for a given \[Theta]", "Section",
 CellChangeTimes->{{3.776508263557572*^9, 
  3.776508290312253*^9}},ExpressionUUID->"2cc81282-9da0-4885-9c78-\
9cf60985c362"],

Cell[CellGroupData[{

Cell[BoxData[
 RowBox[{"Manipulate", "[", 
  RowBox[{
   RowBox[{"Plot", "[", 
    RowBox[{
     RowBox[{
      RowBox[{"Fn", "[", 
       RowBox[{"\[CapitalTheta]", ",", "t"}], "]"}], "//.", 
      "\[VeryThinSpace]", 
      RowBox[{"{", "\[VeryThinSpace]", 
       RowBox[{
        RowBox[{"\[Theta]1", "\[Rule]", 
         RowBox[{"\[Theta]0", "+", " ", 
          FractionBox["A", "\[CapitalLambda]"]}]}], ",", 
        RowBox[{"\[Alpha]", "\[Rule]", 
         FractionBox["1", "13"]}], ",", 
        RowBox[{"\[Beta]", "\[Rule]", "b"}], ",", 
        RowBox[{"A", "\[Rule]", "a"}], ",", 
        RowBox[{"\[CapitalLambda]", "\[Rule]", "\[Lambda]"}], ",", 
        RowBox[{"\[Theta]0", "\[Rule]", "\[CurlyRho]"}]}], "}"}]}], ",", 
     RowBox[{"{", 
      RowBox[{"t", ",", "0", ",", "60"}], "}"}], ",", 
     RowBox[{"PlotRange", "\[Rule]", "All"}], ",", 
     RowBox[{"Filling", "\[Rule]", "Bottom"}]}], "]"}], ",", 
   RowBox[{"{", 
    RowBox[{
     RowBox[{"{", 
      RowBox[{"a", ",", "1", ",", "\"\<A\>\""}], "}"}], ",", "1", ",", "1.5", 
     ",", 
     RowBox[{"Appearance", "\[Rule]", "\"\<Labeled\>\""}]}], "}"}], ",", 
   RowBox[{"{", 
    RowBox[{
     RowBox[{"{", 
      RowBox[{"b", ",", "0.365", ",", "\"\<\[Beta]\>\""}], "}"}], ",", "0.1", 
     ",", "0.6", ",", 
     RowBox[{"Appearance", "\[Rule]", "\"\<Labeled\>\""}]}], "}"}], ",", 
   RowBox[{"{", 
    RowBox[{
     RowBox[{"{", 
      RowBox[{"\[Lambda]", ",", "0.1", ",", "\"\<\[CapitalLambda]\>\""}], 
      "}"}], ",", "0.1", ",", "0.2", ",", 
     RowBox[{"Appearance", "\[Rule]", "\"\<Labeled\>\""}]}], "}"}], ",", 
   RowBox[{"{", 
    RowBox[{
     RowBox[{"{", 
      RowBox[{"\[CurlyRho]", ",", "360", ",", "\"\<\[Theta]0\>\""}], "}"}], 
     ",", "350", ",", "370", ",", 
     RowBox[{"Appearance", "\[Rule]", "\"\<Labeled\>\""}]}], "}"}], ",", 
   RowBox[{"{", 
    RowBox[{
     RowBox[{"{", 
      RowBox[{"\[CapitalTheta]", ",", "380", ",", "\"\<\[Theta]\>\""}], "}"}],
      ",", "\[CurlyRho]", ",", "400", ",", 
     RowBox[{"Appearance", "\[Rule]", "\"\<Labeled\>\""}]}], "}"}], ",", 
   RowBox[{"SaveDefinitions", "\[Rule]", "True"}]}], "]"}]], "Input",
 InitializationCell->True,
 CellChangeTimes->{{3.7764089115682364`*^9, 3.7764090226098704`*^9}, {
   3.776409067331358*^9, 3.7764090779762816`*^9}, {3.7764093224678173`*^9, 
   3.776409572396212*^9}, {3.7764096224396653`*^9, 3.7764096638050823`*^9}, {
   3.7764097399995947`*^9, 3.7764097578963833`*^9}, {3.776409799262787*^9, 
   3.776409833496255*^9}, {3.7764099851317463`*^9, 3.77641000493445*^9}, {
   3.7764100382314663`*^9, 3.7764100440411367`*^9}, {3.7764101549478645`*^9, 
   3.7764101610623755`*^9}, {3.7764104033728766`*^9, 3.776410436862774*^9}, {
   3.7764104795384417`*^9, 3.776410549721387*^9}, 3.776428331968004*^9, {
   3.7764371243388977`*^9, 3.7764371303715067`*^9}, 3.776437164978725*^9, {
   3.77643721536313*^9, 3.7764372240422792`*^9}, 3.7765168427504597`*^9, 
   3.77679868052215*^9, 3.7767987211839066`*^9, {3.780509629098896*^9, 
   3.7805096570845895`*^9}, {3.78055958474354*^9, 
   3.78055960787493*^9}},ExpressionUUID->"6440052f-01b3-4211-94b9-\
5df0a3a02e3b"],

Cell[BoxData[
 TagBox[
  StyleBox[
   DynamicModuleBox[{$CellContext`a$$ = 1, $CellContext`b$$ = 
    0.365, $CellContext`\[CapitalTheta]$$ = 380, $CellContext`\[Lambda]$$ = 
    0.1, $CellContext`\[CurlyRho]$$ = 360, Typeset`show$$ = True, 
    Typeset`bookmarkList$$ = {}, Typeset`bookmarkMode$$ = "Menu", 
    Typeset`animator$$, Typeset`animvar$$ = 1, Typeset`name$$ = 
    "\"untitled\"", Typeset`specs$$ = {{{
       Hold[$CellContext`a$$], 1, "A"}, 1, 1.5}, {{
       Hold[$CellContext`b$$], 0.365, "\[Beta]"}, 0.1, 0.6}, {{
       Hold[$CellContext`\[Lambda]$$], 0.1, "\[CapitalLambda]"}, 0.1, 0.2}, {{
       Hold[$CellContext`\[CurlyRho]$$], 360, "\[Theta]0"}, 350, 370}, {{
       Hold[$CellContext`\[CapitalTheta]$$], 380, "\[Theta]"}, 
      Dynamic[$CellContext`\[CurlyRho]$$], 400}}, Typeset`size$$ = {
    360., {106., 112.}}, Typeset`update$$ = 0, Typeset`initDone$$, 
    Typeset`skipInitDone$$ = False, $CellContext`a$18631$$ = 
    0, $CellContext`b$18632$$ = 0, $CellContext`\[Lambda]$18633$$ = 
    0, $CellContext`\[CurlyRho]$18634$$ = 
    0, $CellContext`\[CapitalTheta]$18635$$ = 0}, 
    DynamicBox[Manipulate`ManipulateBoxes[
     1, StandardForm, 
      "Variables" :> {$CellContext`a$$ = 1, $CellContext`b$$ = 
        0.365, $CellContext`\[CapitalTheta]$$ = 380, $CellContext`\[Lambda]$$ = 
        0.1, $CellContext`\[CurlyRho]$$ = 360}, "ControllerVariables" :> {
        Hold[$CellContext`a$$, $CellContext`a$18631$$, 0], 
        Hold[$CellContext`b$$, $CellContext`b$18632$$, 0], 
        Hold[$CellContext`\[Lambda]$$, $CellContext`\[Lambda]$18633$$, 0], 
        Hold[$CellContext`\[CurlyRho]$$, $CellContext`\[CurlyRho]$18634$$, 0], 
        Hold[$CellContext`\[CapitalTheta]$$, \
$CellContext`\[CapitalTheta]$18635$$, 0]}, 
      "OtherVariables" :> {
       Typeset`show$$, Typeset`bookmarkList$$, Typeset`bookmarkMode$$, 
        Typeset`animator$$, Typeset`animvar$$, Typeset`name$$, 
        Typeset`specs$$, Typeset`size$$, Typeset`update$$, Typeset`initDone$$,
         Typeset`skipInitDone$$}, "Body" :> Plot[
        ReplaceRepeated[
         $CellContext`Fn[$CellContext`\[CapitalTheta]$$, $CellContext`t], \
{$CellContext`\[Theta]1 -> $CellContext`\[Theta]0 + \
$CellContext`A/$CellContext`\[CapitalLambda], $CellContext`\[Alpha] -> 
          1/13, $CellContext`\[Beta] -> $CellContext`b$$, $CellContext`A -> \
$CellContext`a$$, $CellContext`\[CapitalLambda] -> $CellContext`\[Lambda]$$, \
$CellContext`\[Theta]0 -> $CellContext`\[CurlyRho]$$}], {$CellContext`t, 0, 
         60}, PlotRange -> All, Filling -> Bottom], 
      "Specifications" :> {{{$CellContext`a$$, 1, "A"}, 1, 1.5, Appearance -> 
         "Labeled"}, {{$CellContext`b$$, 0.365, "\[Beta]"}, 0.1, 0.6, 
         Appearance -> 
         "Labeled"}, {{$CellContext`\[Lambda]$$, 0.1, "\[CapitalLambda]"}, 
         0.1, 0.2, Appearance -> 
         "Labeled"}, {{$CellContext`\[CurlyRho]$$, 360, "\[Theta]0"}, 350, 
         370, Appearance -> 
         "Labeled"}, {{$CellContext`\[CapitalTheta]$$, 380, "\[Theta]"}, 
         Dynamic[$CellContext`\[CurlyRho]$$], 400, Appearance -> "Labeled"}}, 
      "Options" :> {}, "DefaultOptions" :> {}],
     ImageSizeCache->{411., {209., 215.}},
     SingleEvaluation->True],
    Deinitialization:>None,
    DynamicModuleValues:>{},
    Initialization:>({$CellContext`Fn[
         Pattern[$CellContext`\[Theta], 
          Blank[]], 
         Pattern[$CellContext`t, 
          Blank[]]] := 
       If[$CellContext`\[Theta] > $CellContext`\[Theta]1, 
         If[$CellContext`t > ($CellContext`\[Theta] - \
$CellContext`\[Theta]1)/$CellContext`A, 
           
           Exp[(-$CellContext`\[Alpha]) $CellContext`t - $CellContext`\
\[CapitalLambda] ($CellContext`t - ($CellContext`\[Theta] - $CellContext`\
\[Theta]1)/$CellContext`A) - (($CellContext`\[Beta] \
$CellContext`A)/$CellContext`\[CapitalLambda]) 
            Exp[(-$CellContext`\[CapitalLambda]) ($CellContext`t - \
($CellContext`\[Theta] - $CellContext`\[Theta]1)/$CellContext`A)]], 
           
           Exp[(-($CellContext`\[Alpha] - $CellContext`\[Beta] \
$CellContext`A)) $CellContext`t - $CellContext`\[Beta] ($CellContext`\[Theta] - \
$CellContext`\[Theta]0)]]/((
           E^($CellContext`\[Beta] (-$CellContext`\[Theta] + $CellContext`\
\[Theta]0)) (-1 + 
            E^(((-$CellContext`\[Alpha] + $CellContext`A \
$CellContext`\[Beta]) ($CellContext`\[Theta] - \
$CellContext`\[Theta]1))/$CellContext`A)))/(-$CellContext`\[Alpha] + \
$CellContext`A $CellContext`\[Beta]) + ((
            1/$CellContext`\[CapitalLambda]) $CellContext`G[($CellContext`\
\[Alpha] + $CellContext`\[CapitalLambda])/$CellContext`\[CapitalLambda], \
($CellContext`A $CellContext`\[Beta])/$CellContext`\[CapitalLambda]])/
          E^(($CellContext`\[Alpha] ($CellContext`\[Theta] - $CellContext`\
\[Theta]1))/$CellContext`A)), ($CellContext`\[CapitalLambda] 
          Exp[(-($CellContext`\[Alpha] + $CellContext`\[CapitalLambda])) \
$CellContext`t - ($CellContext`\[Beta] ($CellContext`\[Theta] - $CellContext`\
\[Theta]0)) 
            Exp[(-$CellContext`\[CapitalLambda]) \
$CellContext`t]])/$CellContext`G[($CellContext`\[Alpha] + $CellContext`\
\[CapitalLambda])/$CellContext`\[CapitalLambda], $CellContext`\[Beta] \
($CellContext`\[Theta] - $CellContext`\[Theta]0)]], $CellContext`G[
         Pattern[$CellContext`a, 
          Blank[]], 
         Pattern[$CellContext`\[Xi], 
          Blank[]]] := (Gamma[$CellContext`a] - 
         Gamma[$CellContext`a, \
$CellContext`\[Xi]])/$CellContext`\[Xi]^$CellContext`a}; 
     Typeset`initDone$$ = True),
    SynchronousInitialization->True,
    UndoTrackedVariables:>{Typeset`show$$, Typeset`bookmarkMode$$},
    UnsavedVariables:>{Typeset`initDone$$},
    UntrackedVariables:>{Typeset`size$$}], "Manipulate",
   Deployed->True,
   StripOnInput->False],
  Manipulate`InterpretManipulate[1]]], "Output",
 CellChangeTimes->{3.7765168444704547`*^9, 3.7765182094111147`*^9, 
  3.7765812197623577`*^9, 3.7767066346275516`*^9, 3.776707175031243*^9, 
  3.7767074016663885`*^9, 3.776798725122657*^9, 3.780309166956872*^9, 
  3.7803163720360146`*^9, 3.780422310717002*^9, 3.7804879357233205`*^9, 
  3.7805096604355345`*^9, 3.780514432824829*^9, 3.7805149348928328`*^9, 
  3.78056076143814*^9},ExpressionUUID->"4807522d-91a0-49b9-8501-843b5ba2b458"]
}, Open  ]]
}, Open  ]]
}, Open  ]],

Cell[CellGroupData[{

Cell["Mean age", "Chapter",
 CellChangeTimes->{{3.776411964643236*^9, 3.7764119905884347`*^9}, {
   3.7764284828879137`*^9, 3.776428500361361*^9}, {3.776428914916774*^9, 
   3.7764289253572474`*^9}, {3.776437680426756*^9, 3.776437693866437*^9}, 
   3.7764389439851418`*^9, 3.7764541379835625`*^9, {3.776706828187917*^9, 
   3.776706830194764*^9}},ExpressionUUID->"a78d3a04-98f5-4f39-ae53-\
ce7070c88d9f"],

Cell[CellGroupData[{

Cell["For \[Theta] < \[Theta]1", "Section",
 CellChangeTimes->{{3.776437700204116*^9, 3.776437701066161*^9}, {
   3.776437735580882*^9, 3.7764377556899548`*^9}, 
   3.7764541459180126`*^9},ExpressionUUID->"811db1ff-74a3-4036-9db6-\
315c487ccdfb"],

Cell[CellGroupData[{

Cell[BoxData[
 RowBox[{"B0", " ", "=", " ", 
  RowBox[{
   RowBox[{"Integrate", "[", 
    RowBox[{
     RowBox[{"t", " ", 
      RowBox[{"Exp", "[", 
       RowBox[{
        RowBox[{
         RowBox[{"-", 
          RowBox[{"(", 
           RowBox[{"\[Alpha]", "+", "\[CapitalLambda]"}], ")"}]}], "t"}], " ",
         "-", " ", 
        RowBox[{"\[Beta]", 
         RowBox[{"(", 
          RowBox[{"\[Theta]", "-", "\[Theta]0"}], ")"}], " ", 
         RowBox[{"Exp", "[", 
          RowBox[{
           RowBox[{"-", "\[CapitalLambda]"}], " ", "t"}], "]"}]}]}], "]"}]}], 
     ",", 
     RowBox[{"{", 
      RowBox[{"t", ",", "0", ",", "\[Infinity]"}], "}"}]}], "]"}], "/", 
   "I0"}]}]], "Input",
 CellChangeTimes->{{3.780462741281427*^9, 
  3.7804627417032723`*^9}},ExpressionUUID->"6e1dcdd2-0610-4537-a3e2-\
a00b68e25821"],

Cell[BoxData[
 FractionBox[
  RowBox[{"HypergeometricPFQ", "[", 
   RowBox[{
    RowBox[{"{", 
     RowBox[{
      RowBox[{"1", "+", 
       FractionBox["\[Alpha]", "\[CapitalLambda]"]}], ",", 
      RowBox[{"1", "+", 
       FractionBox["\[Alpha]", "\[CapitalLambda]"]}]}], "}"}], ",", 
    RowBox[{"{", 
     RowBox[{
      RowBox[{"2", "+", 
       FractionBox["\[Alpha]", "\[CapitalLambda]"]}], ",", 
      RowBox[{"2", "+", 
       FractionBox["\[Alpha]", "\[CapitalLambda]"]}]}], "}"}], ",", 
    RowBox[{"\[Beta]", " ", 
     RowBox[{"(", 
      RowBox[{
       RowBox[{"-", "\[Theta]"}], "+", "\[Theta]0"}], ")"}]}]}], "]"}], 
  RowBox[{"I0", " ", 
   SuperscriptBox[
    RowBox[{"(", 
     RowBox[{"\[Alpha]", "+", "\[CapitalLambda]"}], ")"}], "2"]}]]], "Output",
 CellChangeTimes->{3.7804627538301983`*^9, 
  3.780514486439561*^9},ExpressionUUID->"737742ed-bc19-41f5-8f82-\
1fa39ef9f81b"]
}, Open  ]]
}, Open  ]],

Cell[CellGroupData[{

Cell["For \[Theta] > \[Theta]1", "Section",
 CellChangeTimes->{{3.776437875226406*^9, 3.7764379444249496`*^9}, {
  3.776437975121181*^9, 3.776438099883463*^9}, {3.776454156469983*^9, 
  3.7764541666571245`*^9}},ExpressionUUID->"328d4a86-1a09-44c6-a226-\
922c8750058b"],

Cell[CellGroupData[{

Cell["Second part t > (\[Theta]-\[Theta]1)/A", "Subsection",
 CellChangeTimes->{{3.776437875226406*^9, 3.7764379444249496`*^9}, {
  3.776437975121181*^9, 3.776438099883463*^9}, {3.776454156469983*^9, 
  3.77645417691824*^9}, {3.776454856220769*^9, 
  3.7764548619115057`*^9}},ExpressionUUID->"200a26cf-4b03-4423-b796-\
5092098a7d9b"],

Cell[CellGroupData[{

Cell[BoxData[
 RowBox[{"J22", " ", "=", " ", 
  RowBox[{
   RowBox[{"Exp", " ", "[", 
    RowBox[{
     RowBox[{"-", 
      FractionBox["\[Alpha]", "A"]}], 
     RowBox[{"(", 
      RowBox[{"\[Theta]", "-", "\[Theta]1"}], ")"}]}], "]"}], 
   RowBox[{"Assuming", "[", 
    RowBox[{
     RowBox[{"{", 
      RowBox[{
       RowBox[{
        RowBox[{"\[Alpha]", "+", "\[CapitalLambda]"}], ">", "0"}], ",", 
       RowBox[{
        FractionBox[
         RowBox[{"\[Beta]", " ", "A"}], "\[CapitalLambda]"], ">", "0"}], ",", 
       RowBox[{"\[CapitalLambda]", ">", "0"}]}], "}"}], ",", 
     RowBox[{"Integrate", "[", 
      RowBox[{
       RowBox[{"t", " ", 
        RowBox[{"Exp", "[", 
         RowBox[{
          RowBox[{
           RowBox[{"-", 
            RowBox[{"(", 
             RowBox[{"\[Alpha]", "+", "\[CapitalLambda]"}], ")"}]}], " ", 
           "t"}], " ", "-", 
          RowBox[{
           FractionBox[
            RowBox[{"\[Beta]", " ", "A"}], "\[CapitalLambda]"], 
           RowBox[{"Exp", "[", 
            RowBox[{
             RowBox[{"-", "\[CapitalLambda]"}], " ", "t"}], "]"}]}]}], 
         "]"}]}], ",", 
       RowBox[{"{", 
        RowBox[{"t", ",", "0", ",", "\[Infinity]"}], "}"}]}], "]"}]}], 
    "]"}]}]}]], "Input",
 CellChangeTimes->{{3.7767948450335846`*^9, 3.776794878575407*^9}, {
   3.776794933314105*^9, 3.776794938663044*^9}, {3.7767949939524403`*^9, 
   3.776795023893345*^9}, {3.776795070476718*^9, 3.776795082224003*^9}, {
   3.780462828060052*^9, 3.7804628285132046`*^9}, 3.780463459844939*^9, 
   3.780465092321245*^9},ExpressionUUID->"ff942149-91c0-490f-9fb0-\
9bb9cf49422d"],

Cell[BoxData[
 FractionBox[
  RowBox[{
   SuperscriptBox["\[ExponentialE]", 
    RowBox[{"-", 
     FractionBox[
      RowBox[{"\[Alpha]", " ", 
       RowBox[{"(", 
        RowBox[{"\[Theta]", "-", "\[Theta]1"}], ")"}]}], "A"]}]], " ", 
   RowBox[{"HypergeometricPFQ", "[", 
    RowBox[{
     RowBox[{"{", 
      RowBox[{
       RowBox[{"1", "+", 
        FractionBox["\[Alpha]", "\[CapitalLambda]"]}], ",", 
       RowBox[{"1", "+", 
        FractionBox["\[Alpha]", "\[CapitalLambda]"]}]}], "}"}], ",", 
     RowBox[{"{", 
      RowBox[{
       RowBox[{"2", "+", 
        FractionBox["\[Alpha]", "\[CapitalLambda]"]}], ",", 
       RowBox[{"2", "+", 
        FractionBox["\[Alpha]", "\[CapitalLambda]"]}]}], "}"}], ",", 
     RowBox[{"-", 
      FractionBox[
       RowBox[{"A", " ", "\[Beta]"}], "\[CapitalLambda]"]}]}], "]"}]}], 
  SuperscriptBox[
   RowBox[{"(", 
    RowBox[{"\[Alpha]", "+", "\[CapitalLambda]"}], ")"}], "2"]]], "Output",
 CellChangeTimes->{3.7767948934888763`*^9, 3.7767949504063315`*^9, 
  3.7767950303036623`*^9, 3.7767951061523056`*^9, 3.78031707106994*^9, 
  3.7804641765919027`*^9, 3.780465109091796*^9, 
  3.7805144977585*^9},ExpressionUUID->"0442df89-9f9a-4ece-8c6a-24c9ec58d648"]
}, Open  ]]
}, Open  ]],

Cell[CellGroupData[{

Cell["First part t < (\[Theta]-\[Theta]1)/A", "Subsection",
 CellChangeTimes->{{3.7764381974155493`*^9, 3.7764382285348377`*^9}, {
   3.7764383491095667`*^9, 3.77643836104848*^9}, 3.7764541851914964`*^9, {
   3.776454737663742*^9, 3.7764547475700808`*^9}, {3.7764548317607927`*^9, 
   3.77645484913785*^9}},ExpressionUUID->"2ffbfaa1-afd7-438c-be37-\
6ce62fad66f4"],

Cell[CellGroupData[{

Cell[BoxData[
 RowBox[{"J1", " ", "=", 
  RowBox[{
   RowBox[{
    RowBox[{"Exp", "[", 
     RowBox[{
      RowBox[{"-", "\[Beta]"}], 
      RowBox[{"(", 
       RowBox[{"\[Theta]", "-", "\[Theta]0"}], ")"}]}], "]"}], " ", 
    RowBox[{"FullSimplify", "[", 
     RowBox[{"Integrate", "[", 
      RowBox[{
       RowBox[{"t", " ", 
        RowBox[{"Exp", "[", 
         RowBox[{
          RowBox[{"-", 
           RowBox[{"(", 
            RowBox[{"\[Alpha]", " ", "-", " ", 
             RowBox[{"\[Beta]", " ", "A"}]}], ")"}]}], "t"}], "]"}]}], ",", 
       RowBox[{"{", 
        RowBox[{"t", ",", "0", ",", 
         FractionBox[
          RowBox[{"\[Theta]", "-", "\[Theta]1"}], "A"]}], "}"}]}], "]"}], 
     "]"}]}], "/.", 
   RowBox[{
    RowBox[{
     RowBox[{"\[Alpha]", " ", 
      RowBox[{"(", 
       RowBox[{
        RowBox[{"-", "\[Theta]"}], "+", "\[Theta]1"}], ")"}]}], "+", 
     RowBox[{"A", " ", 
      RowBox[{"(", 
       RowBox[{
        RowBox[{"-", "1"}], "+", 
        RowBox[{"\[Beta]", " ", "\[Theta]"}], "-", 
        RowBox[{"\[Beta]", " ", "\[Theta]1"}]}], ")"}]}]}], "\[Rule]", " ", 
    RowBox[{
     RowBox[{
      RowBox[{"(", 
       RowBox[{
        RowBox[{"A", " ", "\[Beta]"}], " ", "-", "\[Alpha]"}], ")"}], 
      RowBox[{"(", 
       RowBox[{"\[Theta]", "-", "\[Theta]1"}], ")"}]}], "-", 
     "A"}]}]}]}]], "Input",
 CellChangeTimes->{{3.7767951543827567`*^9, 3.7767952269702454`*^9}, {
   3.7767952829762344`*^9, 3.7767952856427107`*^9}, 3.780463779841181*^9, {
   3.7804639226420507`*^9, 3.780463963838933*^9}, {3.780463996570265*^9, 
   3.780463997382705*^9}},ExpressionUUID->"e5d69a03-dec6-4fb2-b9f0-\
98abeddef99e"],

Cell[BoxData[
 FractionBox[
  RowBox[{
   SuperscriptBox["\[ExponentialE]", 
    RowBox[{
     RowBox[{"-", "\[Beta]"}], " ", 
     RowBox[{"(", 
      RowBox[{"\[Theta]", "-", "\[Theta]0"}], ")"}]}]], " ", 
   RowBox[{"(", 
    RowBox[{"1", "+", 
     FractionBox[
      RowBox[{
       SuperscriptBox["\[ExponentialE]", 
        FractionBox[
         RowBox[{
          RowBox[{"(", 
           RowBox[{
            RowBox[{"-", "\[Alpha]"}], "+", 
            RowBox[{"A", " ", "\[Beta]"}]}], ")"}], " ", 
          RowBox[{"(", 
           RowBox[{"\[Theta]", "-", "\[Theta]1"}], ")"}]}], "A"]], " ", 
       RowBox[{"(", 
        RowBox[{
         RowBox[{"-", "A"}], "+", 
         RowBox[{
          RowBox[{"(", 
           RowBox[{
            RowBox[{"-", "\[Alpha]"}], "+", 
            RowBox[{"A", " ", "\[Beta]"}]}], ")"}], " ", 
          RowBox[{"(", 
           RowBox[{"\[Theta]", "-", "\[Theta]1"}], ")"}]}]}], ")"}]}], "A"]}],
     ")"}]}], 
  SuperscriptBox[
   RowBox[{"(", 
    RowBox[{"\[Alpha]", "-", 
     RowBox[{"A", " ", "\[Beta]"}]}], ")"}], "2"]]], "Output",
 CellChangeTimes->{{3.776795211621022*^9, 3.7767952287992167`*^9}, 
   3.7767953358989897`*^9, 3.780317050135111*^9, 3.780463797672415*^9, 
   3.7804639657606125`*^9, 3.780464099703747*^9, 3.780483753876848*^9, 
   3.7805144990865507`*^9},ExpressionUUID->"7a5a5d94-1f3f-47ac-9c5b-\
4b7f16cb494e"]
}, Open  ]],

Cell[BoxData[{
 RowBox[{
  RowBox[{
  "The", " ", "parameter", " ", "controling", " ", "the", " ", "shape", " ", 
   "is", " ", "\[Gamma]"}], "=", 
  RowBox[{"\[Beta]", "-", 
   FractionBox["\[Alpha]", "A"]}]}], "\[IndentingNewLine]", 
 RowBox[{"J1", " ", "=", " ", 
  RowBox[{
   FractionBox[
    SuperscriptBox["\[ExponentialE]", 
     RowBox[{
      RowBox[{"-", "\[Beta]"}], " ", 
      RowBox[{"(", 
       RowBox[{"\[Theta]", "-", "\[Theta]0"}], ")"}]}]], 
    SuperscriptBox[
     RowBox[{"(", 
      RowBox[{"\[Alpha]", "-", 
       RowBox[{"A", " ", "\[Beta]"}]}], ")"}], "2"]], " ", 
   RowBox[{"(", 
    RowBox[{"1", "+", 
     RowBox[{
      SuperscriptBox["\[ExponentialE]", 
       RowBox[{"\[Gamma]", " ", 
        RowBox[{"(", 
         RowBox[{"\[Theta]", "-", "\[Theta]1"}], ")"}]}]], " ", 
      RowBox[{"(", 
       RowBox[{
        RowBox[{"-", "1"}], "+", 
        RowBox[{"\[Gamma]", " ", 
         RowBox[{"(", 
          RowBox[{"\[Theta]", "-", "\[Theta]1"}], ")"}]}]}], ")"}]}]}], 
    ")"}]}]}]}], "Text",
 CellChangeTimes->{{3.780484807420829*^9, 3.780484868588704*^9}, {
  3.7804861657127433`*^9, 3.7804861793341236`*^9}, {3.7804863725779643`*^9, 
  3.780486467465478*^9}},ExpressionUUID->"c000d6a2-7ed1-4c4f-9beb-\
b6b9af979965"]
}, Open  ]],

Cell[CellGroupData[{

Cell["Together", "Subsection",
 CellChangeTimes->{{3.7804666907541285`*^9, 3.780466702957122*^9}, {
  3.7805596397330713`*^9, 
  3.780559647516816*^9}},ExpressionUUID->"f06ec256-c77d-4014-927e-\
58faa4334509"],

Cell[CellGroupData[{

Cell["Results are first repeated for easy initialization", "Subsubsection",
 CellChangeTimes->{{3.7805596510165687`*^9, 
  3.7805596926568947`*^9}},ExpressionUUID->"57b2f3e2-f8c1-4c57-bb83-\
83e683d5c1ca"],

Cell[BoxData[{
 RowBox[{
  RowBox[{
   RowBox[{"BO", "=", " ", 
    FractionBox[
     RowBox[{"HypergeometricPFQ", "[", 
      RowBox[{
       RowBox[{"{", 
        RowBox[{
         RowBox[{"1", "+", 
          FractionBox["\[Alpha]", "\[CapitalLambda]"]}], ",", 
         RowBox[{"1", "+", 
          FractionBox["\[Alpha]", "\[CapitalLambda]"]}]}], "}"}], ",", 
       RowBox[{"{", 
        RowBox[{
         RowBox[{"2", "+", 
          FractionBox["\[Alpha]", "\[CapitalLambda]"]}], ",", 
         RowBox[{"2", "+", 
          FractionBox["\[Alpha]", "\[CapitalLambda]"]}]}], "}"}], ",", 
       RowBox[{"\[Beta]", " ", 
        RowBox[{"(", 
         RowBox[{
          RowBox[{"-", "\[Theta]"}], "+", "\[Theta]0"}], ")"}]}]}], "]"}], 
     RowBox[{"I0", " ", 
      SuperscriptBox[
       RowBox[{"(", 
        RowBox[{"\[Alpha]", "+", "\[CapitalLambda]"}], ")"}], "2"]}]]}], 
   ";"}], "\[IndentingNewLine]"}], "\[IndentingNewLine]", 
 RowBox[{
  RowBox[{"J22", " ", "=", " ", 
   FractionBox[
    RowBox[{
     SuperscriptBox["\[ExponentialE]", 
      RowBox[{"-", 
       FractionBox[
        RowBox[{"\[Alpha]", " ", 
         RowBox[{"(", 
          RowBox[{"\[Theta]", "-", "\[Theta]1"}], ")"}]}], "A"]}]], " ", 
     RowBox[{"HypergeometricPFQ", "[", 
      RowBox[{
       RowBox[{"{", 
        RowBox[{
         RowBox[{"1", "+", 
          FractionBox["\[Alpha]", "\[CapitalLambda]"]}], ",", 
         RowBox[{"1", "+", 
          FractionBox["\[Alpha]", "\[CapitalLambda]"]}]}], "}"}], ",", 
       RowBox[{"{", 
        RowBox[{
         RowBox[{"2", "+", 
          FractionBox["\[Alpha]", "\[CapitalLambda]"]}], ",", 
         RowBox[{"2", "+", 
          FractionBox["\[Alpha]", "\[CapitalLambda]"]}]}], "}"}], ",", 
       RowBox[{"-", 
        FractionBox[
         RowBox[{"A", " ", "\[Beta]"}], "\[CapitalLambda]"]}]}], "]"}]}], 
    SuperscriptBox[
     RowBox[{"(", 
      RowBox[{"\[Alpha]", "+", "\[CapitalLambda]"}], ")"}], "2"]]}], 
  ";"}], "\[IndentingNewLine]", 
 RowBox[{
  RowBox[{"J1", " ", "=", " ", 
   FractionBox[
    RowBox[{
     SuperscriptBox["\[ExponentialE]", 
      RowBox[{
       RowBox[{"-", "\[Beta]"}], " ", 
       RowBox[{"(", 
        RowBox[{"\[Theta]", "-", "\[Theta]0"}], ")"}]}]], " ", 
     RowBox[{"(", 
      RowBox[{"1", "+", 
       FractionBox[
        RowBox[{
         SuperscriptBox["\[ExponentialE]", 
          FractionBox[
           RowBox[{
            RowBox[{"(", 
             RowBox[{
              RowBox[{"-", "\[Alpha]"}], "+", 
              RowBox[{"A", " ", "\[Beta]"}]}], ")"}], " ", 
            RowBox[{"(", 
             RowBox[{"\[Theta]", "-", "\[Theta]1"}], ")"}]}], "A"]], " ", 
         RowBox[{"(", 
          RowBox[{
           RowBox[{"-", "A"}], "+", 
           RowBox[{
            RowBox[{"(", 
             RowBox[{
              RowBox[{"-", "\[Alpha]"}], "+", 
              RowBox[{"A", " ", "\[Beta]"}]}], ")"}], " ", 
            RowBox[{"(", 
             RowBox[{"\[Theta]", "-", "\[Theta]1"}], ")"}]}]}], ")"}]}], 
        "A"]}], ")"}]}], 
    SuperscriptBox[
     RowBox[{"(", 
      RowBox[{"\[Alpha]", "-", 
       RowBox[{"A", " ", "\[Beta]"}]}], ")"}], "2"]]}], ";"}]}], "Input",
 InitializationCell->True,
 CellChangeTimes->{{3.7805596807215133`*^9, 
  3.7805597661540203`*^9}},ExpressionUUID->"9d0b4dfa-f0ea-4362-ace5-\
51e3fe15721d"]
}, Closed]],

Cell[CellGroupData[{

Cell["Definition of the mean age", "Subsubsection",
 CellChangeTimes->{{3.7805605379338264`*^9, 
  3.7805605653927717`*^9}},ExpressionUUID->"754de79c-2935-4a6c-bb47-\
3f088351eb09"],

Cell[BoxData[
 RowBox[{
  RowBox[{"Mage", "[", "\[Theta]_", "]"}], ":=", 
  RowBox[{"If", "[", 
   RowBox[{
    RowBox[{"\[Theta]", "<", " ", "\[Theta]1"}], ",", "B0", ",", 
    FractionBox[
     RowBox[{"J1", " ", "+", "  ", "J22", " ", "+", " ", 
      RowBox[{
       FractionBox[
        RowBox[{"\[Theta]", "-", "\[Theta]1"}], "A"], "I2"}]}], 
     RowBox[{"I1", " ", "+", "  ", "I2"}]]}], "]"}]}]], "Input",
 CellChangeTimes->{{3.776707124891887*^9, 3.776707150865041*^9}, {
   3.776707181086768*^9, 3.7767071815834846`*^9}, 3.7767072644111376`*^9, {
   3.776707459653244*^9, 3.7767074991525383`*^9}, {3.7767955498626795`*^9, 
   3.776795595063859*^9}, {3.776795865887041*^9, 3.7767960789062743`*^9}, {
   3.7767961572584915`*^9, 3.776796186591724*^9}, {3.7804635812658377`*^9, 
   3.7804635817189336`*^9}, {3.7804636156386538`*^9, 
   3.7804636164042053`*^9}, {3.7804636769331665`*^9, 3.780463677495562*^9}, {
   3.780464007585353*^9, 3.7804640479463625`*^9}, {3.780465063368255*^9, 
   3.780465068586624*^9}, 
   3.7804667367937107`*^9},ExpressionUUID->"66d5669c-5161-4653-919f-\
eee583243123"]
}, Open  ]],

Cell[CellGroupData[{

Cell["Mean age as a function of \[Theta]", "Subsubsection",
 CellChangeTimes->{{3.780560578224886*^9, 
  3.7805606091416135`*^9}},ExpressionUUID->"9286bf45-5db9-4cd0-91dd-\
8c6e1c3f06bd"],

Cell[BoxData[{
 RowBox[{"The", " ", "mean", " ", "age", " ", "curve", " ", "is", " ", 
  "plotted", " ", "together", " ", "with", " ", "the", " ", 
  FractionBox[
   RowBox[{"\[Theta]", "-", "\[Theta]1"}], "A"], 
  "curve"}], "\[IndentingNewLine]", 
 RowBox[{"Tech", ":", " ", 
  RowBox[{
  "Definitions", " ", "need", " ", "to", " ", "be", " ", "repeated", " ", 
   "within", " ", "the", " ", "cell", " ", "for", " ", "the", " ", "CDF", " ",
    "player", " ", "to", " ", "know", " ", 
   RowBox[{"them", "."}]}]}]}], "Text",
 CellChangeTimes->{{3.7805622243872423`*^9, 3.7805622606666217`*^9}, {
  3.7805631004080496`*^9, 
  3.7805631794373856`*^9}},ExpressionUUID->"45b6e362-0c36-463b-959a-\
bf15cb35087a"],

Cell[CellGroupData[{

Cell[BoxData[{
 RowBox[{
  RowBox[{
   RowBox[{"G", "[", 
    RowBox[{"a_", ",", "\[Xi]_"}], "]"}], ":=", " ", 
   RowBox[{
    SuperscriptBox["\[Xi]", 
     RowBox[{"-", "a"}]], " ", 
    RowBox[{"(", 
     RowBox[{
      RowBox[{"Gamma", "[", "a", "]"}], "-", 
      RowBox[{"Gamma", "[", 
       RowBox[{"a", ",", "\[Xi]"}], "]"}]}], ")"}]}]}], 
  ";"}], "\[IndentingNewLine]", 
 RowBox[{
  RowBox[{
   RowBox[{"I0", " ", "=", " ", 
    RowBox[{
     FractionBox["1", "\[CapitalLambda]"], 
     RowBox[{"G", "[", 
      RowBox[{
       FractionBox[
        RowBox[{"\[Alpha]", "+", "\[CapitalLambda]"}], "\[CapitalLambda]"], 
       ",", 
       RowBox[{"\[Beta]", " ", 
        RowBox[{"(", 
         RowBox[{"\[Theta]", "-", "\[Theta]0"}], ")"}]}]}], "]"}]}]}], ";"}], 
  "\[IndentingNewLine]"}], "\[IndentingNewLine]", 
 RowBox[{
  RowBox[{"I2", " ", "=", " ", 
   RowBox[{
    SuperscriptBox["\[ExponentialE]", 
     RowBox[{"-", 
      FractionBox[
       RowBox[{"\[Alpha]", " ", 
        RowBox[{"(", 
         RowBox[{"\[Theta]", "-", "\[Theta]1"}], ")"}]}], "A"]}]], 
    FractionBox["1", "\[CapitalLambda]"], 
    RowBox[{"G", "[", 
     RowBox[{
      FractionBox[
       RowBox[{"\[Alpha]", "+", "\[CapitalLambda]"}], "\[CapitalLambda]"], 
      ",", 
      FractionBox[
       RowBox[{"A", " ", "\[Beta]"}], "\[CapitalLambda]"]}], "]"}]}]}], 
  ";"}], "\[IndentingNewLine]", 
 RowBox[{
  RowBox[{
   RowBox[{"I1", " ", "=", " ", 
    FractionBox[
     RowBox[{
      SuperscriptBox["\[ExponentialE]", 
       RowBox[{"\[Beta]", " ", 
        RowBox[{"(", 
         RowBox[{
          RowBox[{"-", "\[Theta]"}], "+", "\[Theta]0"}], ")"}]}]], " ", 
      RowBox[{"(", 
       RowBox[{
        RowBox[{"-", "1"}], "+", 
        SuperscriptBox["\[ExponentialE]", 
         FractionBox[
          RowBox[{
           RowBox[{"(", 
            RowBox[{
             RowBox[{"-", "\[Alpha]"}], "+", 
             RowBox[{"A", " ", "\[Beta]"}]}], ")"}], " ", 
           RowBox[{"(", 
            RowBox[{"\[Theta]", "-", "\[Theta]1"}], ")"}]}], "A"]]}], ")"}]}], 
     RowBox[{
      RowBox[{"-", "\[Alpha]"}], "+", 
      RowBox[{"A", " ", "\[Beta]"}]}]]}], ";"}], 
  "\[IndentingNewLine]"}], "\[IndentingNewLine]", 
 RowBox[{
  RowBox[{"B0", "=", 
   FractionBox[
    RowBox[{"HypergeometricPFQ", "[", 
     RowBox[{
      RowBox[{"{", 
       RowBox[{
        RowBox[{"1", "+", 
         FractionBox["\[Alpha]", "\[CapitalLambda]"]}], ",", 
        RowBox[{"1", "+", 
         FractionBox["\[Alpha]", "\[CapitalLambda]"]}]}], "}"}], ",", 
      RowBox[{"{", 
       RowBox[{
        RowBox[{"2", "+", 
         FractionBox["\[Alpha]", "\[CapitalLambda]"]}], ",", 
        RowBox[{"2", "+", 
         FractionBox["\[Alpha]", "\[CapitalLambda]"]}]}], "}"}], ",", 
      RowBox[{"\[Beta]", " ", 
       RowBox[{"(", 
        RowBox[{
         RowBox[{"-", "\[Theta]"}], "+", "\[Theta]0"}], ")"}]}]}], "]"}], 
    RowBox[{"I0", " ", 
     SuperscriptBox[
      RowBox[{"(", 
       RowBox[{"\[Alpha]", "+", "\[CapitalLambda]"}], ")"}], "2"]}]]}], 
  ";"}], "\[IndentingNewLine]", 
 RowBox[{
  RowBox[{"J1", " ", "=", " ", 
   FractionBox[
    RowBox[{
     SuperscriptBox["\[ExponentialE]", 
      RowBox[{
       RowBox[{"-", "\[Beta]"}], " ", 
       RowBox[{"(", 
        RowBox[{"\[Theta]", "-", "\[Theta]0"}], ")"}]}]], " ", 
     RowBox[{"(", 
      RowBox[{"1", "+", 
       FractionBox[
        RowBox[{
         SuperscriptBox["\[ExponentialE]", 
          FractionBox[
           RowBox[{
            RowBox[{"(", 
             RowBox[{
              RowBox[{"-", "\[Alpha]"}], "+", 
              RowBox[{"A", " ", "\[Beta]"}]}], ")"}], " ", 
            RowBox[{"(", 
             RowBox[{"\[Theta]", "-", "\[Theta]1"}], ")"}]}], "A"]], " ", 
         RowBox[{"(", 
          RowBox[{
           RowBox[{"-", "A"}], "+", 
           RowBox[{
            RowBox[{"(", 
             RowBox[{
              RowBox[{"-", "\[Alpha]"}], "+", 
              RowBox[{"A", " ", "\[Beta]"}]}], ")"}], " ", 
            RowBox[{"(", 
             RowBox[{"\[Theta]", "-", "\[Theta]1"}], ")"}]}]}], ")"}]}], 
        "A"]}], ")"}]}], 
    SuperscriptBox[
     RowBox[{"(", 
      RowBox[{"\[Alpha]", "-", 
       RowBox[{"A", " ", "\[Beta]"}]}], ")"}], "2"]]}], 
  ";"}], "\[IndentingNewLine]", 
 RowBox[{
  RowBox[{
   RowBox[{"J22", " ", "=", " ", 
    FractionBox[
     RowBox[{
      SuperscriptBox["\[ExponentialE]", 
       RowBox[{"-", 
        FractionBox[
         RowBox[{"\[Alpha]", " ", 
          RowBox[{"(", 
           RowBox[{"\[Theta]", "-", "\[Theta]1"}], ")"}]}], "A"]}]], " ", 
      RowBox[{"HypergeometricPFQ", "[", 
       RowBox[{
        RowBox[{"{", 
         RowBox[{
          RowBox[{"1", "+", 
           FractionBox["\[Alpha]", "\[CapitalLambda]"]}], ",", 
          RowBox[{"1", "+", 
           FractionBox["\[Alpha]", "\[CapitalLambda]"]}]}], "}"}], ",", 
        RowBox[{"{", 
         RowBox[{
          RowBox[{"2", "+", 
           FractionBox["\[Alpha]", "\[CapitalLambda]"]}], ",", 
          RowBox[{"2", "+", 
           FractionBox["\[Alpha]", "\[CapitalLambda]"]}]}], "}"}], ",", 
        RowBox[{"-", 
         FractionBox[
          RowBox[{"A", " ", "\[Beta]"}], "\[CapitalLambda]"]}]}], "]"}]}], 
     SuperscriptBox[
      RowBox[{"(", 
       RowBox[{"\[Alpha]", "+", "\[CapitalLambda]"}], ")"}], "2"]]}], ";"}], 
  "\[IndentingNewLine]"}], "\[IndentingNewLine]", 
 RowBox[{
  RowBox[{
   RowBox[{"Mage", "[", "\[Theta]_", "]"}], ":=", 
   RowBox[{"If", "[", 
    RowBox[{
     RowBox[{"\[Theta]", "<", " ", "\[Theta]1"}], ",", "B0", ",", 
     FractionBox[
      RowBox[{"J1", " ", "+", "  ", "J22", " ", "+", " ", 
       RowBox[{
        FractionBox[
         RowBox[{"\[Theta]", "-", "\[Theta]1"}], "A"], "I2"}]}], 
      RowBox[{"I1", " ", "+", "  ", "I2"}]]}], "]"}]}], 
  "\[IndentingNewLine]"}], "\[IndentingNewLine]", 
 RowBox[{"Manipulate", "[", 
  RowBox[{
   RowBox[{"Plot", "[", 
    RowBox[{
     RowBox[{
      RowBox[{"{", 
       RowBox[{
        RowBox[{"Mage", "[", "\[Theta]", "]"}], " ", ",", 
        FractionBox[
         RowBox[{"\[Theta]", "-", "\[Theta]1"}], "A"]}], "}"}], "//.", 
      "\[VeryThinSpace]", 
      RowBox[{"{", "\[VeryThinSpace]", 
       RowBox[{
        RowBox[{"\[Theta]1", "\[Rule]", 
         RowBox[{"\[Theta]0", "+", " ", 
          FractionBox["A", "\[CapitalLambda]"]}]}], ",", 
        RowBox[{"\[Alpha]", "\[Rule]", 
         FractionBox["1", "13"]}], ",", 
        RowBox[{"\[Beta]", "\[Rule]", "b"}], ",", 
        RowBox[{"A", "\[Rule]", "a"}], ",", 
        RowBox[{"\[CapitalLambda]", "\[Rule]", "\[Lambda]"}], ",", 
        RowBox[{"\[Theta]0", "\[Rule]", "\[CurlyRho]"}]}], "}"}]}], ",", 
     RowBox[{"{", 
      RowBox[{"\[Theta]", ",", "\[CurlyRho]", ",", "440"}], "}"}], ",", 
     RowBox[{"PlotRange", "\[Rule]", "All"}]}], "]"}], ",", 
   RowBox[{"{", 
    RowBox[{
     RowBox[{"{", 
      RowBox[{"a", ",", "1", ",", "\"\<A\>\""}], "}"}], ",", "0.5", ",", 
     "1.5", ",", 
     RowBox[{"Appearance", "\[Rule]", "\"\<Labeled\>\""}]}], "}"}], ",", 
   RowBox[{"{", 
    RowBox[{
     RowBox[{"{", 
      RowBox[{"b", ",", "0.365", ",", "\"\<\[Beta]\>\""}], "}"}], ",", "0.1", 
     ",", "0.6", ",", 
     RowBox[{"Appearance", "\[Rule]", "\"\<Labeled\>\""}]}], "}"}], ",", 
   RowBox[{"{", 
    RowBox[{
     RowBox[{"{", 
      RowBox[{"\[Lambda]", ",", "0.1", ",", "\"\<\[CapitalLambda]\>\""}], 
      "}"}], ",", "0.1", ",", "0.3", ",", 
     RowBox[{"Appearance", "\[Rule]", "\"\<Labeled\>\""}]}], "}"}], ",", 
   RowBox[{"{", 
    RowBox[{
     RowBox[{"{", 
      RowBox[{"\[CurlyRho]", ",", "360", ",", "\"\<\[Theta]0\>\""}], "}"}], 
     ",", "350", ",", "370", ",", 
     RowBox[{"Appearance", "\[Rule]", "\"\<Labeled\>\""}]}], "}"}], ",", 
   RowBox[{"SaveDefinitions", "\[Rule]", "True"}]}], "]"}]}], "Input",
 InitializationCell->True,
 CellChangeTimes->{{3.7767075459138055`*^9, 3.776707675218894*^9}, {
   3.7767077065539846`*^9, 3.776707745141926*^9}, {3.7767984387513466`*^9, 
   3.7767985013495703`*^9}, {3.776798914981083*^9, 3.7767989397868867`*^9}, {
   3.7767994793604794`*^9, 3.7767994846024637`*^9}, {3.780379076902196*^9, 
   3.7803790777163115`*^9}, 3.7804645419791803`*^9, {3.780562141073993*^9, 
   3.78056215216856*^9}, {3.780562967437479*^9, 
   3.780562980882306*^9}},ExpressionUUID->"25b6c3b0-f5a8-423d-99a8-\
cd5c7159ddb9"],

Cell[BoxData[
 TagBox[
  StyleBox[
   DynamicModuleBox[{$CellContext`a$$ = 1, $CellContext`b$$ = 
    0.365, $CellContext`\[Lambda]$$ = 0.1, $CellContext`\[CurlyRho]$$ = 360, 
    Typeset`show$$ = True, Typeset`bookmarkList$$ = {}, 
    Typeset`bookmarkMode$$ = "Menu", Typeset`animator$$, Typeset`animvar$$ = 
    1, Typeset`name$$ = "\"untitled\"", Typeset`specs$$ = {{{
       Hold[$CellContext`a$$], 1, "A"}, 0.5, 1.5}, {{
       Hold[$CellContext`b$$], 0.365, "\[Beta]"}, 0.1, 0.6}, {{
       Hold[$CellContext`\[Lambda]$$], 0.1, "\[CapitalLambda]"}, 0.1, 0.3}, {{
       Hold[$CellContext`\[CurlyRho]$$], 360, "\[Theta]0"}, 350, 370}}, 
    Typeset`size$$ = {360., {106., 112.}}, Typeset`update$$ = 0, 
    Typeset`initDone$$, Typeset`skipInitDone$$ = 
    False, $CellContext`a$39158$$ = 0, $CellContext`b$39159$$ = 
    0, $CellContext`\[Lambda]$39160$$ = 0, $CellContext`\[CurlyRho]$39161$$ = 
    0}, 
    DynamicBox[Manipulate`ManipulateBoxes[
     1, StandardForm, 
      "Variables" :> {$CellContext`a$$ = 1, $CellContext`b$$ = 
        0.365, $CellContext`\[Lambda]$$ = 0.1, $CellContext`\[CurlyRho]$$ = 
        360}, "ControllerVariables" :> {
        Hold[$CellContext`a$$, $CellContext`a$39158$$, 0], 
        Hold[$CellContext`b$$, $CellContext`b$39159$$, 0], 
        Hold[$CellContext`\[Lambda]$$, $CellContext`\[Lambda]$39160$$, 0], 
        Hold[$CellContext`\[CurlyRho]$$, $CellContext`\[CurlyRho]$39161$$, 
         0]}, "OtherVariables" :> {
       Typeset`show$$, Typeset`bookmarkList$$, Typeset`bookmarkMode$$, 
        Typeset`animator$$, Typeset`animvar$$, Typeset`name$$, 
        Typeset`specs$$, Typeset`size$$, Typeset`update$$, Typeset`initDone$$,
         Typeset`skipInitDone$$}, "Body" :> Plot[
        ReplaceRepeated[{
          $CellContext`Mage[$CellContext`\[Theta]], ($CellContext`\[Theta] - \
$CellContext`\[Theta]1)/$CellContext`A}, {$CellContext`\[Theta]1 -> \
$CellContext`\[Theta]0 + $CellContext`A/$CellContext`\[CapitalLambda], \
$CellContext`\[Alpha] -> 
          1/13, $CellContext`\[Beta] -> $CellContext`b$$, $CellContext`A -> \
$CellContext`a$$, $CellContext`\[CapitalLambda] -> $CellContext`\[Lambda]$$, \
$CellContext`\[Theta]0 -> $CellContext`\[CurlyRho]$$}], \
{$CellContext`\[Theta], $CellContext`\[CurlyRho]$$, 440}, PlotRange -> All], 
      "Specifications" :> {{{$CellContext`a$$, 1, "A"}, 0.5, 1.5, Appearance -> 
         "Labeled"}, {{$CellContext`b$$, 0.365, "\[Beta]"}, 0.1, 0.6, 
         Appearance -> 
         "Labeled"}, {{$CellContext`\[Lambda]$$, 0.1, "\[CapitalLambda]"}, 
         0.1, 0.3, Appearance -> 
         "Labeled"}, {{$CellContext`\[CurlyRho]$$, 360, "\[Theta]0"}, 350, 
         370, Appearance -> "Labeled"}}, "Options" :> {}, 
      "DefaultOptions" :> {}],
     ImageSizeCache->{411., {195., 201.}},
     SingleEvaluation->True],
    Deinitialization:>None,
    DynamicModuleValues:>{},
    Initialization:>({$CellContext`Mage[
         Pattern[$CellContext`\[Theta], 
          Blank[]]] := 
       If[$CellContext`\[Theta] < $CellContext`\[Theta]1, $CellContext`B0, \
($CellContext`J1 + $CellContext`J22 + (($CellContext`\[Theta] - $CellContext`\
\[Theta]1)/$CellContext`A) $CellContext`I2)/($CellContext`I1 + \
$CellContext`I2)], $CellContext`B0 = ((($CellContext`\[Beta] ($CellContext`\
\[Theta] - $CellContext`\[Theta]0))^(($CellContext`\[Alpha] + $CellContext`\
\[CapitalLambda])/$CellContext`\[CapitalLambda]) \
$CellContext`\[CapitalLambda]) 
         HypergeometricPFQ[{
           1 + $CellContext`\[Alpha]/$CellContext`\[CapitalLambda], 
            1 + $CellContext`\[Alpha]/$CellContext`\[CapitalLambda]}, {
           2 + $CellContext`\[Alpha]/$CellContext`\[CapitalLambda], 
            2 + $CellContext`\[Alpha]/$CellContext`\[CapitalLambda]}, \
$CellContext`\[Beta] (-$CellContext`\[Theta] + \
$CellContext`\[Theta]0)])/(($CellContext`\[Alpha] + $CellContext`\
\[CapitalLambda])^2 (
         Gamma[($CellContext`\[Alpha] + \
$CellContext`\[CapitalLambda])/$CellContext`\[CapitalLambda]] - 
         Gamma[($CellContext`\[Alpha] + \
$CellContext`\[CapitalLambda])/$CellContext`\[CapitalLambda], $CellContext`\
\[Beta] ($CellContext`\[Theta] - $CellContext`\[Theta]0)])), $CellContext`J1 = \
(1 + (E^(((-$CellContext`\[Alpha] + $CellContext`A $CellContext`\[Beta]) \
($CellContext`\[Theta] - $CellContext`\[Theta]1))/$CellContext`A) \
(-$CellContext`A + (-$CellContext`\[Alpha] + $CellContext`A $CellContext`\
\[Beta]) ($CellContext`\[Theta] - $CellContext`\[Theta]1)))/$CellContext`A)/(
        E^($CellContext`\[Beta] ($CellContext`\[Theta] - \
$CellContext`\[Theta]0)) ($CellContext`\[Alpha] - $CellContext`A \
$CellContext`\[Beta])^2), $CellContext`J22 = 
       HypergeometricPFQ[{
          1 + $CellContext`\[Alpha]/$CellContext`\[CapitalLambda], 
           1 + $CellContext`\[Alpha]/$CellContext`\[CapitalLambda]}, {
          2 + $CellContext`\[Alpha]/$CellContext`\[CapitalLambda], 
           2 + $CellContext`\[Alpha]/$CellContext`\[CapitalLambda]}, \
-(($CellContext`A $CellContext`\[Beta])/$CellContext`\[CapitalLambda])]/(
        E^(($CellContext`\[Alpha] ($CellContext`\[Theta] - $CellContext`\
\[Theta]1))/$CellContext`A) ($CellContext`\[Alpha] + $CellContext`\
\[CapitalLambda])^2), $CellContext`I2 = (
         Gamma[($CellContext`\[Alpha] + \
$CellContext`\[CapitalLambda])/$CellContext`\[CapitalLambda]] - 
         Gamma[($CellContext`\[Alpha] + \
$CellContext`\[CapitalLambda])/$CellContext`\[CapitalLambda], ($CellContext`A \
$CellContext`\[Beta])/$CellContext`\[CapitalLambda]])/((
         E^(($CellContext`\[Alpha] ($CellContext`\[Theta] - $CellContext`\
\[Theta]1))/$CellContext`A) (($CellContext`A \
$CellContext`\[Beta])/$CellContext`\[CapitalLambda])^(($CellContext`\[Alpha] + \
$CellContext`\[CapitalLambda])/$CellContext`\[CapitalLambda])) $CellContext`\
\[CapitalLambda]), $CellContext`I1 = (
         E^($CellContext`\[Beta] (-$CellContext`\[Theta] + $CellContext`\
\[Theta]0)) (-1 + 
          E^(((-$CellContext`\[Alpha] + $CellContext`A $CellContext`\[Beta]) \
($CellContext`\[Theta] - \
$CellContext`\[Theta]1))/$CellContext`A)))/(-$CellContext`\[Alpha] + \
$CellContext`A $CellContext`\[Beta])}; Typeset`initDone$$ = True),
    SynchronousInitialization->True,
    UndoTrackedVariables:>{Typeset`show$$, Typeset`bookmarkMode$$},
    UnsavedVariables:>{Typeset`initDone$$},
    UntrackedVariables:>{Typeset`size$$}], "Manipulate",
   Deployed->True,
   StripOnInput->False],
  Manipulate`InterpretManipulate[1]]], "Output",
 CellChangeTimes->{{3.7767077310170193`*^9, 3.776707746719026*^9}, 
   3.776796107413*^9, 3.776798506551612*^9, 3.7767989411830854`*^9, 
   3.7767995953491597`*^9, 3.780316743569078*^9, 3.7803170983710327`*^9, 
   3.7803790815597286`*^9, 3.780464120624076*^9, 3.7804641877803335`*^9, 
   3.780464859523575*^9, 3.780465128705978*^9, 3.7804667918420825`*^9, 
   3.7805144472975473`*^9, 3.7805145130230455`*^9, 3.78056076203185*^9, 
   3.7805629992455273`*^9},ExpressionUUID->"010f8b2f-f465-4c3c-9110-\
56bdb7d4d6b9"]
}, Open  ]],

Cell[TextData[{
 "There are essentially two cases for  \[Theta] > \[Theta]1, see below, one \
for which  J1 dominates the numerator and I1 dominates the denominator, for \
small values of \[Beta], and one for which ",
 Cell[BoxData[
  RowBox[{"J22", " ", "+", " ", 
   RowBox[{
    FractionBox[
     RowBox[{"\[Theta]", "-", "\[Theta]1"}], "A"], "I2"}]}]],
  InitializationCell->True,
  CellChangeTimes->{{3.776707124891887*^9, 3.776707150865041*^9}, {
    3.776707181086768*^9, 3.7767071815834846`*^9}, 3.7767072644111376`*^9, {
    3.776707459653244*^9, 3.7767074991525383`*^9}, {3.7767955498626795`*^9, 
    3.776795595063859*^9}, {3.776795865887041*^9, 3.7767960789062743`*^9}, {
    3.7767961572584915`*^9, 3.776796186591724*^9}, {3.7804635812658377`*^9, 
    3.7804635817189336`*^9}, {3.7804636156386538`*^9, 
    3.7804636164042053`*^9}, {3.7804636769331665`*^9, 3.780463677495562*^9}, {
    3.780464007585353*^9, 3.7804640479463625`*^9}, {3.780465063368255*^9, 
    3.780465068586624*^9}, 3.7804667367937107`*^9},ExpressionUUID->
  "75a21b3a-9b50-4c27-9522-b8e3711eb34a"],
 " dominates the numerator and I2 dominates the denominator. In the first \
case, the mean age is   ",
 Cell[BoxData[
  FormBox[
   FractionBox["J1", "I1"], TraditionalForm]],ExpressionUUID->
  "f6196a91-e7b6-43aa-90e5-1b142045ce7d"],
 Cell[BoxData[
  RowBox[{"=", 
   RowBox[{
    FractionBox["1", "A"], 
    RowBox[{"(", 
     RowBox[{
      RowBox[{"-", 
       FractionBox["1", "\[Gamma]"]}], "+", "\[Theta]", "-", "\[Theta]1", "+", 
      FractionBox[
       RowBox[{"\[Theta]", "-", "\[Theta]1"}], 
       RowBox[{
        RowBox[{"-", "1"}], "+", 
        SuperscriptBox["\[ExponentialE]", 
         RowBox[{"\[Gamma]", " ", 
          RowBox[{"(", 
           RowBox[{"\[Theta]", "-", "\[Theta]1"}], ")"}]}]]}]]}]}]}]}]],
  CellChangeTimes->{3.780488891905812*^9, 3.7804889254335594`*^9, 
   3.7804889628760586`*^9, 3.78048914867153*^9},ExpressionUUID->
  "c8c9a72f-c69a-4d3a-af67-36ed8315fb4f"],
 "). If \[Gamma]>0, the age behaves at large \[Theta] as ",
 Cell[BoxData[
  FractionBox["1", "A"]],
  CellChangeTimes->{3.780488891905812*^9, 3.7804889254335594`*^9, 
   3.7804889628760586`*^9, 3.78048914867153*^9},ExpressionUUID->
  "3e2e1937-8026-4be1-86f4-57b65bf2e14f"],
 Cell[BoxData[
  SuperscriptBox["\[ExponentialE]", 
   RowBox[{
    RowBox[{"-", "\[Gamma]"}], " ", 
    RowBox[{"(", 
     RowBox[{"\[Theta]", "-", "\[Theta]1"}], ")"}]}]]],
  CellChangeTimes->{3.780488891905812*^9, 3.7804889254335594`*^9, 
   3.7804889628760586`*^9, 3.78048914867153*^9},ExpressionUUID->
  "3aadb3e8-5078-4121-9683-ec87183d400c"],
 "(\[Theta]-\[Theta]1), and if \[Gamma]<0, the age behaves as ",
 Cell[BoxData[
  FractionBox["1", "A"]],
  CellChangeTimes->{3.780488891905812*^9, 3.7804889254335594`*^9, 
   3.7804889628760586`*^9, 3.78048914867153*^9},ExpressionUUID->
  "216697ab-942b-4d0e-92e5-adf9e2799b4c"],
 " (\[Theta]-\[Theta]1). In the second case, the mean age is ",
 Cell[BoxData[
  FractionBox[
   RowBox[{"J22", " ", "+", " ", 
    RowBox[{
     FractionBox[
      RowBox[{"\[Theta]", "-", "\[Theta]1"}], "A"], "I2"}]}], "I2"]],
  InitializationCell->True,
  CellChangeTimes->{{3.776707124891887*^9, 3.776707150865041*^9}, {
    3.776707181086768*^9, 3.7767071815834846`*^9}, 3.7767072644111376`*^9, {
    3.776707459653244*^9, 3.7767074991525383`*^9}, {3.7767955498626795`*^9, 
    3.776795595063859*^9}, {3.776795865887041*^9, 3.7767960789062743`*^9}, {
    3.7767961572584915`*^9, 3.776796186591724*^9}, {3.7804635812658377`*^9, 
    3.7804635817189336`*^9}, {3.7804636156386538`*^9, 
    3.7804636164042053`*^9}, {3.7804636769331665`*^9, 3.780463677495562*^9}, {
    3.780464007585353*^9, 3.7804640479463625`*^9}, {3.780465063368255*^9, 
    3.780465068586624*^9}, 3.7804667367937107`*^9},ExpressionUUID->
  "5ca8ff3d-5b57-47c8-b66a-0301e94b42ee"],
 "= ",
 Cell[BoxData[
  FractionBox[
   RowBox[{"\[Theta]", "-", "\[Theta]1"}], "A"]],
  InitializationCell->True,
  CellChangeTimes->{{3.776707124891887*^9, 3.776707150865041*^9}, {
    3.776707181086768*^9, 3.7767071815834846`*^9}, 3.7767072644111376`*^9, {
    3.776707459653244*^9, 3.7767074991525383`*^9}, {3.7767955498626795`*^9, 
    3.776795595063859*^9}, {3.776795865887041*^9, 3.7767960789062743`*^9}, {
    3.7767961572584915`*^9, 3.776796186591724*^9}, {3.7804635812658377`*^9, 
    3.7804635817189336`*^9}, {3.7804636156386538`*^9, 
    3.7804636164042053`*^9}, {3.7804636769331665`*^9, 3.780463677495562*^9}, {
    3.780464007585353*^9, 3.7804640479463625`*^9}, {3.780465063368255*^9, 
    3.780465068586624*^9}, 3.7804667367937107`*^9},ExpressionUUID->
  "49685dfe-5166-46a6-bb7d-70161717fc54"],
 "+ ",
 Cell[BoxData[
  FormBox[
   FractionBox["J22", "I2"], TraditionalForm]],ExpressionUUID->
  "c396a09b-0a4e-4c52-91c1-58716dc1b33e"],
 " where ",
 Cell[BoxData[
  FormBox[
   FractionBox["J22", "I2"], TraditionalForm]],ExpressionUUID->
  "8ca467ac-eb30-4711-ae9e-22a9ff489bf3"],
 " does not depend on \[Theta].\nTherefore, except in the small \[Beta] and \
positive \[Gamma] case, the  mean age growth with \[Theta] follows the mean \
diabatic ascent rate A."
}], "Text",
 CellChangeTimes->{{3.7767986644243603`*^9, 3.7767986650809746`*^9}, {
   3.7767988402148185`*^9, 3.7767988759723625`*^9}, {3.776799408973694*^9, 
   3.7767994643190603`*^9}, {3.77679967157759*^9, 3.7767996800957193`*^9}, {
   3.776799723992629*^9, 3.776799736229636*^9}, 3.7803796762925053`*^9, 
   3.780465312861783*^9, {3.78048850293266*^9, 3.780488540909853*^9}, 
   3.7804886296129184`*^9, {3.780488671207015*^9, 3.780488812940175*^9}, {
   3.7804891772950974`*^9, 3.7804892175749645`*^9}, {3.780489248824792*^9, 
   3.780489260302182*^9}, {3.780489345068508*^9, 3.7804893798825345`*^9}, {
   3.7804894199165287`*^9, 3.7804899151556497`*^9}, {3.7804900533762827`*^9, 
   3.780490054761485*^9}},ExpressionUUID->"ad386b5b-589f-4589-9096-\
21a4a2415beb"]
}, Open  ]],

Cell[CellGroupData[{

Cell["Analysis of the contributions to the mean age", "Subsubsection",
 CellChangeTimes->{{3.780484469605501*^9, 
  3.7804845051900697`*^9}},ExpressionUUID->"fd34bb06-f856-4e0a-ad7c-\
7c5ab2a0670c"],

Cell["Terms in the denominator", "Item",
 CellChangeTimes->{{3.7804845207061596`*^9, 3.7804845477606087`*^9}, 
   3.7804845788307695`*^9},ExpressionUUID->"57189437-f234-459d-b259-\
65d1ac87b850"],

Cell[CellGroupData[{

Cell[BoxData[
 RowBox[{"Manipulate", "[", 
  RowBox[{
   RowBox[{"Plot", "[", 
    RowBox[{
     RowBox[{"{", 
      RowBox[{
       RowBox[{"I1", " ", "//.", 
        RowBox[{"{", "\[VeryThinSpace]", 
         RowBox[{
          RowBox[{"\[Theta]1", "\[Rule]", 
           RowBox[{"\[Theta]0", "+", " ", 
            FractionBox["A", "\[CapitalLambda]"]}]}], ",", 
          RowBox[{"\[Alpha]", "\[Rule]", 
           FractionBox["1", "13"]}], ",", 
          RowBox[{"\[Beta]", "\[Rule]", "b"}], ",", 
          RowBox[{"A", "\[Rule]", "a"}], ",", 
          RowBox[{"\[CapitalLambda]", "\[Rule]", "\[Lambda]"}], ",", 
          RowBox[{"\[Theta]0", "\[Rule]", "\[CurlyRho]"}]}], "}"}]}], " ", 
       ",", 
       RowBox[{"I2", "//.", 
        RowBox[{"{", "\[VeryThinSpace]", 
         RowBox[{
          RowBox[{"\[Theta]1", "\[Rule]", 
           RowBox[{"\[Theta]0", "+", " ", 
            FractionBox["A", "\[CapitalLambda]"]}]}], ",", 
          RowBox[{"\[Alpha]", "\[Rule]", 
           FractionBox["1", "13"]}], ",", 
          RowBox[{"\[Beta]", "\[Rule]", "b"}], ",", 
          RowBox[{"A", "\[Rule]", "a"}], ",", 
          RowBox[{"\[CapitalLambda]", "\[Rule]", "\[Lambda]"}], ",", 
          RowBox[{"\[Theta]0", "\[Rule]", "\[CurlyRho]"}]}], "}"}]}]}], "}"}],
      "\[VeryThinSpace]", ",", 
     RowBox[{"{", 
      RowBox[{"\[Theta]", ",", 
       RowBox[{"\[CurlyRho]", "+", 
        FractionBox["a", "\[Lambda]"]}], ",", "440"}], "}"}], ",", 
     RowBox[{"PlotStyle", "\[Rule]", 
      RowBox[{"{", 
       RowBox[{"Blue", ",", "Red"}], "}"}]}], ",", 
     RowBox[{"PlotRange", "\[Rule]", "All"}]}], "]"}], ",", 
   RowBox[{"{", 
    RowBox[{
     RowBox[{"{", 
      RowBox[{"a", ",", "1", ",", "\"\<A\>\""}], "}"}], ",", "0.5", ",", 
     "1.5", ",", 
     RowBox[{"Appearance", "\[Rule]", "\"\<Labeled\>\""}]}], "}"}], ",", 
   RowBox[{"{", 
    RowBox[{
     RowBox[{"{", 
      RowBox[{"b", ",", "0.325", ",", "\"\<\[Beta]\>\""}], "}"}], ",", "0.1", 
     ",", "0.6", ",", 
     RowBox[{"Appearance", "\[Rule]", "\"\<Labeled\>\""}]}], "}"}], ",", 
   RowBox[{"{", 
    RowBox[{
     RowBox[{"{", 
      RowBox[{"\[Lambda]", ",", "0.1", ",", "\"\<\[CapitalLambda]\>\""}], 
      "}"}], ",", "0.1", ",", "0.3", ",", 
     RowBox[{"Appearance", "\[Rule]", "\"\<Labeled\>\""}]}], "}"}], ",", 
   RowBox[{"{", 
    RowBox[{
     RowBox[{"{", 
      RowBox[{"\[CurlyRho]", ",", "360", ",", "\"\<\[Theta]0\>\""}], "}"}], 
     ",", "350", ",", "370", ",", 
     RowBox[{"Appearance", "\[Rule]", "\"\<Labeled\>\""}]}], "}"}], ",", 
   RowBox[{"SaveDefinitions", "\[Rule]", "True"}]}], "]"}]], "Input",
 InitializationCell->True,
 CellChangeTimes->{{3.7804798991360397`*^9, 3.7804799063059225`*^9}, {
   3.7804805211075406`*^9, 3.780480534586799*^9}, 3.780480672571566*^9, {
   3.7804808132411313`*^9, 3.780480813638914*^9}, {3.7804811823449626`*^9, 
   3.780481273301751*^9}, {3.7804814255113325`*^9, 3.7804814610599194`*^9}, {
   3.7804815418910084`*^9, 3.780481550671954*^9}, {3.7804815828464794`*^9, 
   3.780481626658333*^9}, {3.78048173565948*^9, 3.780481746748113*^9}, {
   3.780481785245996*^9, 3.780481787657611*^9}, {3.7804821165820208`*^9, 
   3.780482192964161*^9}, {3.780482597322678*^9, 3.780482658310687*^9}, 
   3.780482697946928*^9, {3.780482775514388*^9, 3.7804828388570065`*^9}, {
   3.7804829217414017`*^9, 3.780482923999106*^9}, {3.780483955759941*^9, 
   3.7804839722094927`*^9}, {3.780559836747164*^9, 
   3.7805598457267637`*^9}},ExpressionUUID->"71e3c7df-a753-4e76-8c67-\
87df5939c0a5"],

Cell[BoxData[
 TagBox[
  StyleBox[
   DynamicModuleBox[{$CellContext`a$$ = 1, $CellContext`b$$ = 
    0.325, $CellContext`\[Lambda]$$ = 0.1, $CellContext`\[CurlyRho]$$ = 360, 
    Typeset`show$$ = True, Typeset`bookmarkList$$ = {}, 
    Typeset`bookmarkMode$$ = "Menu", Typeset`animator$$, Typeset`animvar$$ = 
    1, Typeset`name$$ = "\"untitled\"", Typeset`specs$$ = {{{
       Hold[$CellContext`a$$], 1, "A"}, 0.5, 1.5}, {{
       Hold[$CellContext`b$$], 0.325, "\[Beta]"}, 0.1, 0.6}, {{
       Hold[$CellContext`\[Lambda]$$], 0.1, "\[CapitalLambda]"}, 0.1, 0.3}, {{
       Hold[$CellContext`\[CurlyRho]$$], 360, "\[Theta]0"}, 350, 370}}, 
    Typeset`size$$ = {360., {108., 113.}}, Typeset`update$$ = 0, 
    Typeset`initDone$$, Typeset`skipInitDone$$ = 
    False, $CellContext`a$19211$$ = 0, $CellContext`b$19212$$ = 
    0, $CellContext`\[Lambda]$19213$$ = 0, $CellContext`\[CurlyRho]$19214$$ = 
    0}, 
    DynamicBox[Manipulate`ManipulateBoxes[
     1, StandardForm, 
      "Variables" :> {$CellContext`a$$ = 1, $CellContext`b$$ = 
        0.325, $CellContext`\[Lambda]$$ = 0.1, $CellContext`\[CurlyRho]$$ = 
        360}, "ControllerVariables" :> {
        Hold[$CellContext`a$$, $CellContext`a$19211$$, 0], 
        Hold[$CellContext`b$$, $CellContext`b$19212$$, 0], 
        Hold[$CellContext`\[Lambda]$$, $CellContext`\[Lambda]$19213$$, 0], 
        Hold[$CellContext`\[CurlyRho]$$, $CellContext`\[CurlyRho]$19214$$, 
         0]}, "OtherVariables" :> {
       Typeset`show$$, Typeset`bookmarkList$$, Typeset`bookmarkMode$$, 
        Typeset`animator$$, Typeset`animvar$$, Typeset`name$$, 
        Typeset`specs$$, Typeset`size$$, Typeset`update$$, Typeset`initDone$$,
         Typeset`skipInitDone$$}, "Body" :> Plot[{
         ReplaceRepeated[$CellContext`I1, {$CellContext`\[Theta]1 -> \
$CellContext`\[Theta]0 + $CellContext`A/$CellContext`\[CapitalLambda], \
$CellContext`\[Alpha] -> 
           1/13, $CellContext`\[Beta] -> $CellContext`b$$, $CellContext`A -> \
$CellContext`a$$, $CellContext`\[CapitalLambda] -> $CellContext`\[Lambda]$$, \
$CellContext`\[Theta]0 -> $CellContext`\[CurlyRho]$$}], 
         ReplaceRepeated[$CellContext`I2, {$CellContext`\[Theta]1 -> \
$CellContext`\[Theta]0 + $CellContext`A/$CellContext`\[CapitalLambda], \
$CellContext`\[Alpha] -> 
           1/13, $CellContext`\[Beta] -> $CellContext`b$$, $CellContext`A -> \
$CellContext`a$$, $CellContext`\[CapitalLambda] -> $CellContext`\[Lambda]$$, \
$CellContext`\[Theta]0 -> $CellContext`\[CurlyRho]$$}]}, {$CellContext`\
\[Theta], $CellContext`\[CurlyRho]$$ + \
$CellContext`a$$/$CellContext`\[Lambda]$$, 440}, PlotStyle -> {Blue, Red}, 
        PlotRange -> All], 
      "Specifications" :> {{{$CellContext`a$$, 1, "A"}, 0.5, 1.5, Appearance -> 
         "Labeled"}, {{$CellContext`b$$, 0.325, "\[Beta]"}, 0.1, 0.6, 
         Appearance -> 
         "Labeled"}, {{$CellContext`\[Lambda]$$, 0.1, "\[CapitalLambda]"}, 
         0.1, 0.3, Appearance -> 
         "Labeled"}, {{$CellContext`\[CurlyRho]$$, 360, "\[Theta]0"}, 350, 
         370, Appearance -> "Labeled"}}, "Options" :> {}, 
      "DefaultOptions" :> {}],
     ImageSizeCache->{411., {197., 203.}},
     SingleEvaluation->True],
    Deinitialization:>None,
    DynamicModuleValues:>{},
    Initialization:>({$CellContext`I1 = (
         E^($CellContext`\[Beta] (-$CellContext`\[Theta] + $CellContext`\
\[Theta]0)) (-1 + 
          E^(((-$CellContext`\[Alpha] + $CellContext`A $CellContext`\[Beta]) \
($CellContext`\[Theta] - \
$CellContext`\[Theta]1))/$CellContext`A)))/(-$CellContext`\[Alpha] + \
$CellContext`A $CellContext`\[Beta]), $CellContext`I2 = (
         Gamma[($CellContext`\[Alpha] + \
$CellContext`\[CapitalLambda])/$CellContext`\[CapitalLambda]] - 
         Gamma[($CellContext`\[Alpha] + \
$CellContext`\[CapitalLambda])/$CellContext`\[CapitalLambda], ($CellContext`A \
$CellContext`\[Beta])/$CellContext`\[CapitalLambda]])/(
        E^(($CellContext`\[Alpha] ($CellContext`\[Theta] - $CellContext`\
\[Theta]1))/$CellContext`A) (($CellContext`A \
$CellContext`\[Beta])/$CellContext`\[CapitalLambda])^(($CellContext`\[Alpha] + \
$CellContext`\[CapitalLambda])/$CellContext`\[CapitalLambda]) $CellContext`\
\[CapitalLambda])}; Typeset`initDone$$ = True),
    SynchronousInitialization->True,
    UndoTrackedVariables:>{Typeset`show$$, Typeset`bookmarkMode$$},
    UnsavedVariables:>{Typeset`initDone$$},
    UntrackedVariables:>{Typeset`size$$}], "Manipulate",
   Deployed->True,
   StripOnInput->False],
  Manipulate`InterpretManipulate[1]]], "Output",
 CellChangeTimes->{3.780484638235656*^9, 3.780559848796311*^9, 
  3.7805607625019035`*^9},ExpressionUUID->"b8d0cbdf-e06c-497f-be15-\
7f9e677c0f0c"]
}, Open  ]],

Cell["\<\
The red curve (I2) dominates the blue curve (I1) except at low  values of \
\[Beta]\
\>", "Text",
 CellChangeTimes->{{3.780484648539747*^9, 3.7804846938857*^9}, {
  3.7804886075995398`*^9, 
  3.7804886233195267`*^9}},ExpressionUUID->"55a6f005-a642-497f-8d03-\
44d5e83c1130"],

Cell["Terms in the numerator", "Item",
 CellChangeTimes->{{3.7804846142164564`*^9, 
  3.7804846280645156`*^9}},ExpressionUUID->"b4c4df69-cec5-403f-a98b-\
a63bc8b4add9"],

Cell[CellGroupData[{

Cell[BoxData[
 RowBox[{"Manipulate", "[", 
  RowBox[{
   RowBox[{"Plot", "[", 
    RowBox[{
     RowBox[{"{", 
      RowBox[{
       RowBox[{"J1", " ", "//.", 
        RowBox[{"{", "\[VeryThinSpace]", 
         RowBox[{
          RowBox[{"\[Theta]1", "\[Rule]", 
           RowBox[{"\[Theta]0", "+", " ", 
            FractionBox["A", "\[CapitalLambda]"]}]}], ",", 
          RowBox[{"\[Alpha]", "\[Rule]", 
           FractionBox["1", "13"]}], ",", 
          RowBox[{"\[Beta]", "\[Rule]", "b"}], ",", 
          RowBox[{"A", "\[Rule]", "a"}], ",", 
          RowBox[{"\[CapitalLambda]", "\[Rule]", "\[Lambda]"}], ",", 
          RowBox[{"\[Theta]0", "\[Rule]", "\[CurlyRho]"}]}], "}"}]}], " ", 
       ",", 
       RowBox[{"J22", "//.", " ", 
        RowBox[{"{", "\[VeryThinSpace]", 
         RowBox[{
          RowBox[{"\[Theta]1", "\[Rule]", 
           RowBox[{"\[Theta]0", "+", " ", 
            FractionBox["A", "\[CapitalLambda]"]}]}], ",", 
          RowBox[{"\[Alpha]", "\[Rule]", 
           FractionBox["1", "13"]}], ",", 
          RowBox[{"\[Beta]", "\[Rule]", "b"}], ",", 
          RowBox[{"A", "\[Rule]", "a"}], ",", 
          RowBox[{"\[CapitalLambda]", "\[Rule]", "\[Lambda]"}], ",", 
          RowBox[{"\[Theta]0", "\[Rule]", "\[CurlyRho]"}]}], "}"}]}], ",", 
       RowBox[{
        RowBox[{
         FractionBox[
          RowBox[{"\[Theta]", "-", "\[Theta]1"}], "A"], "I2"}], "//.", " ", 
        RowBox[{"{", "\[VeryThinSpace]", 
         RowBox[{
          RowBox[{"\[Theta]1", "\[Rule]", 
           RowBox[{"\[Theta]0", "+", " ", 
            FractionBox["A", "\[CapitalLambda]"]}]}], ",", 
          RowBox[{"\[Alpha]", "\[Rule]", 
           FractionBox["1", "13"]}], ",", 
          RowBox[{"\[Beta]", "\[Rule]", "b"}], ",", 
          RowBox[{"A", "\[Rule]", "a"}], ",", 
          RowBox[{"\[CapitalLambda]", "\[Rule]", "\[Lambda]"}], ",", 
          RowBox[{"\[Theta]0", "\[Rule]", "\[CurlyRho]"}]}], "}"}]}]}], "}"}],
      "\[VeryThinSpace]", ",", 
     RowBox[{"{", 
      RowBox[{"\[Theta]", ",", 
       RowBox[{"\[CurlyRho]", "+", 
        FractionBox["a", "\[Lambda]"]}], ",", "440"}], "}"}], ",", 
     RowBox[{"PlotStyle", "\[Rule]", 
      RowBox[{"{", 
       RowBox[{"Blue", ",", "Red", ",", "Green"}], "}"}]}], ",", 
     RowBox[{"PlotRange", "\[Rule]", "All"}]}], "]"}], ",", 
   RowBox[{"{", 
    RowBox[{
     RowBox[{"{", 
      RowBox[{"a", ",", "1", ",", "\"\<A\>\""}], "}"}], ",", "0.5", ",", 
     "1.5", ",", 
     RowBox[{"Appearance", "\[Rule]", "\"\<Labeled\>\""}]}], "}"}], ",", 
   RowBox[{"{", 
    RowBox[{
     RowBox[{"{", 
      RowBox[{"b", ",", "0.325", ",", "\"\<\[Beta]\>\""}], "}"}], ",", "0.1", 
     ",", "0.6", ",", 
     RowBox[{"Appearance", "\[Rule]", "\"\<Labeled\>\""}]}], "}"}], ",", 
   RowBox[{"{", 
    RowBox[{
     RowBox[{"{", 
      RowBox[{"\[Lambda]", ",", "0.1", ",", "\"\<\[CapitalLambda]\>\""}], 
      "}"}], ",", "0.1", ",", "0.3", ",", 
     RowBox[{"Appearance", "\[Rule]", "\"\<Labeled\>\""}]}], "}"}], ",", 
   RowBox[{"{", 
    RowBox[{
     RowBox[{"{", 
      RowBox[{"\[CurlyRho]", ",", "360", ",", "\"\<\[Theta]0\>\""}], "}"}], 
     ",", "350", ",", "370", ",", 
     RowBox[{"Appearance", "\[Rule]", "\"\<Labeled\>\""}]}], "}"}], ",", 
   RowBox[{"SaveDefinitions", "\[Rule]", "True"}]}], "]"}]], "Input",
 InitializationCell->True,
 CellChangeTimes->{{3.7804798991360397`*^9, 3.7804799063059225`*^9}, {
   3.7804805211075406`*^9, 3.780480534586799*^9}, 3.780480672571566*^9, {
   3.7804808132411313`*^9, 3.780480813638914*^9}, {3.7804811823449626`*^9, 
   3.780481273301751*^9}, {3.7804814255113325`*^9, 3.7804814610599194`*^9}, {
   3.7804815418910084`*^9, 3.780481550671954*^9}, {3.7804815828464794`*^9, 
   3.780481626658333*^9}, {3.78048173565948*^9, 3.780481746748113*^9}, {
   3.780481785245996*^9, 3.780481787657611*^9}, {3.7804821165820208`*^9, 
   3.780482192964161*^9}, {3.780482597322678*^9, 3.780482658310687*^9}, 
   3.780482697946928*^9, {3.780482775514388*^9, 3.7804828388570065`*^9}, {
   3.7804829217414017`*^9, 3.780482923999106*^9}, {3.7804836635209436`*^9, 
   3.780483666114169*^9}, {3.7804837112992973`*^9, 3.7804837225368443`*^9}, {
   3.7804837862432623`*^9, 3.780483841181738*^9}, {3.780559874160319*^9, 
   3.7805598854622965`*^9}},ExpressionUUID->"855011be-1371-410a-824d-\
e1f78fdcdf25"],

Cell[BoxData[
 TagBox[
  StyleBox[
   DynamicModuleBox[{$CellContext`a$$ = 1, $CellContext`b$$ = 
    0.325, $CellContext`\[Lambda]$$ = 0.1, $CellContext`\[CurlyRho]$$ = 360, 
    Typeset`show$$ = True, Typeset`bookmarkList$$ = {}, 
    Typeset`bookmarkMode$$ = "Menu", Typeset`animator$$, Typeset`animvar$$ = 
    1, Typeset`name$$ = "\"untitled\"", Typeset`specs$$ = {{{
       Hold[$CellContext`a$$], 1, "A"}, 0.5, 1.5}, {{
       Hold[$CellContext`b$$], 0.325, "\[Beta]"}, 0.1, 0.6}, {{
       Hold[$CellContext`\[Lambda]$$], 0.1, "\[CapitalLambda]"}, 0.1, 0.3}, {{
       Hold[$CellContext`\[CurlyRho]$$], 360, "\[Theta]0"}, 350, 370}}, 
    Typeset`size$$ = {360., {109., 113.}}, Typeset`update$$ = 0, 
    Typeset`initDone$$, Typeset`skipInitDone$$ = 
    False, $CellContext`a$19238$$ = 0, $CellContext`b$19239$$ = 
    0, $CellContext`\[Lambda]$19240$$ = 0, $CellContext`\[CurlyRho]$19241$$ = 
    0}, 
    DynamicBox[Manipulate`ManipulateBoxes[
     1, StandardForm, 
      "Variables" :> {$CellContext`a$$ = 1, $CellContext`b$$ = 
        0.325, $CellContext`\[Lambda]$$ = 0.1, $CellContext`\[CurlyRho]$$ = 
        360}, "ControllerVariables" :> {
        Hold[$CellContext`a$$, $CellContext`a$19238$$, 0], 
        Hold[$CellContext`b$$, $CellContext`b$19239$$, 0], 
        Hold[$CellContext`\[Lambda]$$, $CellContext`\[Lambda]$19240$$, 0], 
        Hold[$CellContext`\[CurlyRho]$$, $CellContext`\[CurlyRho]$19241$$, 
         0]}, "OtherVariables" :> {
       Typeset`show$$, Typeset`bookmarkList$$, Typeset`bookmarkMode$$, 
        Typeset`animator$$, Typeset`animvar$$, Typeset`name$$, 
        Typeset`specs$$, Typeset`size$$, Typeset`update$$, Typeset`initDone$$,
         Typeset`skipInitDone$$}, "Body" :> Plot[{
         ReplaceRepeated[$CellContext`J1, {$CellContext`\[Theta]1 -> \
$CellContext`\[Theta]0 + $CellContext`A/$CellContext`\[CapitalLambda], \
$CellContext`\[Alpha] -> 
           1/13, $CellContext`\[Beta] -> $CellContext`b$$, $CellContext`A -> \
$CellContext`a$$, $CellContext`\[CapitalLambda] -> $CellContext`\[Lambda]$$, \
$CellContext`\[Theta]0 -> $CellContext`\[CurlyRho]$$}], 
         ReplaceRepeated[$CellContext`J22, {$CellContext`\[Theta]1 -> \
$CellContext`\[Theta]0 + $CellContext`A/$CellContext`\[CapitalLambda], \
$CellContext`\[Alpha] -> 
           1/13, $CellContext`\[Beta] -> $CellContext`b$$, $CellContext`A -> \
$CellContext`a$$, $CellContext`\[CapitalLambda] -> $CellContext`\[Lambda]$$, \
$CellContext`\[Theta]0 -> $CellContext`\[CurlyRho]$$}], 
         ReplaceRepeated[(($CellContext`\[Theta] - \
$CellContext`\[Theta]1)/$CellContext`A) $CellContext`I2, {$CellContext`\
\[Theta]1 -> $CellContext`\[Theta]0 + $CellContext`A/$CellContext`\
\[CapitalLambda], $CellContext`\[Alpha] -> 
           1/13, $CellContext`\[Beta] -> $CellContext`b$$, $CellContext`A -> \
$CellContext`a$$, $CellContext`\[CapitalLambda] -> $CellContext`\[Lambda]$$, \
$CellContext`\[Theta]0 -> $CellContext`\[CurlyRho]$$}]}, {$CellContext`\
\[Theta], $CellContext`\[CurlyRho]$$ + \
$CellContext`a$$/$CellContext`\[Lambda]$$, 440}, 
        PlotStyle -> {Blue, Red, Green}, PlotRange -> All], 
      "Specifications" :> {{{$CellContext`a$$, 1, "A"}, 0.5, 1.5, Appearance -> 
         "Labeled"}, {{$CellContext`b$$, 0.325, "\[Beta]"}, 0.1, 0.6, 
         Appearance -> 
         "Labeled"}, {{$CellContext`\[Lambda]$$, 0.1, "\[CapitalLambda]"}, 
         0.1, 0.3, Appearance -> 
         "Labeled"}, {{$CellContext`\[CurlyRho]$$, 360, "\[Theta]0"}, 350, 
         370, Appearance -> "Labeled"}}, "Options" :> {}, 
      "DefaultOptions" :> {}],
     ImageSizeCache->{411., {197., 203.}},
     SingleEvaluation->True],
    Deinitialization:>None,
    DynamicModuleValues:>{},
    Initialization:>({$CellContext`J1 = (
         1 + (E^(((-$CellContext`\[Alpha] + $CellContext`A $CellContext`\
\[Beta]) ($CellContext`\[Theta] - $CellContext`\[Theta]1))/$CellContext`A) \
(-$CellContext`A + (-$CellContext`\[Alpha] + $CellContext`A $CellContext`\
\[Beta]) ($CellContext`\[Theta] - $CellContext`\[Theta]1)))/$CellContext`A)/(
        E^($CellContext`\[Beta] ($CellContext`\[Theta] - \
$CellContext`\[Theta]0)) ($CellContext`\[Alpha] - $CellContext`A \
$CellContext`\[Beta])^2), $CellContext`J22 = 
       HypergeometricPFQ[{
          1 + $CellContext`\[Alpha]/$CellContext`\[CapitalLambda], 
           1 + $CellContext`\[Alpha]/$CellContext`\[CapitalLambda]}, {
          2 + $CellContext`\[Alpha]/$CellContext`\[CapitalLambda], 
           2 + $CellContext`\[Alpha]/$CellContext`\[CapitalLambda]}, \
-(($CellContext`A $CellContext`\[Beta])/$CellContext`\[CapitalLambda])]/(
        E^(($CellContext`\[Alpha] ($CellContext`\[Theta] - $CellContext`\
\[Theta]1))/$CellContext`A) ($CellContext`\[Alpha] + $CellContext`\
\[CapitalLambda])^2), $CellContext`I2 = (
         Gamma[($CellContext`\[Alpha] + \
$CellContext`\[CapitalLambda])/$CellContext`\[CapitalLambda]] - 
         Gamma[($CellContext`\[Alpha] + \
$CellContext`\[CapitalLambda])/$CellContext`\[CapitalLambda], ($CellContext`A \
$CellContext`\[Beta])/$CellContext`\[CapitalLambda]])/(
        E^(($CellContext`\[Alpha] ($CellContext`\[Theta] - $CellContext`\
\[Theta]1))/$CellContext`A) (($CellContext`A \
$CellContext`\[Beta])/$CellContext`\[CapitalLambda])^(($CellContext`\[Alpha] + \
$CellContext`\[CapitalLambda])/$CellContext`\[CapitalLambda]) $CellContext`\
\[CapitalLambda])}; Typeset`initDone$$ = True),
    SynchronousInitialization->True,
    UndoTrackedVariables:>{Typeset`show$$, Typeset`bookmarkMode$$},
    UnsavedVariables:>{Typeset`initDone$$},
    UntrackedVariables:>{Typeset`size$$}], "Manipulate",
   Deployed->True,
   StripOnInput->False],
  Manipulate`InterpretManipulate[1]]], "Output",
 CellChangeTimes->{3.780484642930981*^9, 3.7804868267665267`*^9, 
  3.780560762580043*^9},ExpressionUUID->"5b837fa4-f0f7-4d55-9b59-\
316816067546"]
}, Open  ]],

Cell["\<\
Thechnical note : repeating the rule was found necessary for Plot to account \
for the PlotStyle\
\>", "Text",
 CellChangeTimes->{{3.7804899241784678`*^9, 3.780490002520483*^9}, 
   3.78055989876606*^9},ExpressionUUID->"59be3c0d-9b66-4d37-8f37-\
f15bc9488127"]
}, Closed]]
}, Open  ]]
}, Open  ]]
}, Open  ]],

Cell[CellGroupData[{

Cell["Modal age", "Chapter",
 CellChangeTimes->{{3.776411964643236*^9, 3.7764119905884347`*^9}, {
   3.7764284828879137`*^9, 3.776428500361361*^9}, {3.776428914916774*^9, 
   3.7764289253572474`*^9}, {3.776437680426756*^9, 3.776437693866437*^9}, 
   3.7764389439851418`*^9, 3.7764541379835625`*^9, {3.776706828187917*^9, 
   3.776706830194764*^9}, {3.7767088195687776`*^9, 
   3.7767088206201754`*^9}},ExpressionUUID->"b65f26ef-55e2-4a84-ac82-\
b23e4dc95dde"],

Cell["\<\
The cases for the modal age are described above in the discussion of the age \
section at a given \[Theta] . Rendered here using Max and Min.\
\>", "Text",
 CellChangeTimes->{{3.7767998105941267`*^9, 3.7767998271726685`*^9}, 
   3.7804900671813717`*^9, {3.7804902230121393`*^9, 3.780490246837474*^9}, {
   3.7804902770471087`*^9, 3.7804903275421124`*^9}, {3.7805090292640915`*^9, 
   3.7805090571152844`*^9}},ExpressionUUID->"b06bcadf-1c7e-49f6-b81a-\
6b0e8929ebb1"],

Cell[BoxData[
 RowBox[{
  RowBox[{"ModalAge", "[", "\[Theta]_", "]"}], " ", ":=", 
  RowBox[{
   RowBox[{"Max", "[", 
    RowBox[{
     FractionBox[
      RowBox[{"\[Theta]", "-", "\[Theta]1"}], "A"], ",", "0"}], "]"}], "+", 
   RowBox[{"Max", "[", 
    RowBox[{
     RowBox[{
      FractionBox["1", "\[CapitalLambda]"], 
      RowBox[{"Log", "[", 
       FractionBox[
        RowBox[{"\[Beta]", " ", 
         RowBox[{"Min", "[", 
          RowBox[{"A", ",", 
           RowBox[{"\[CapitalLambda]", 
            RowBox[{"(", 
             RowBox[{"\[Theta]", "-", "\[Theta]0"}], ")"}]}]}], "]"}]}], 
        RowBox[{"\[Alpha]", "+", "\[CapitalLambda]"}]], "]"}]}], ",", "0"}], 
    "]"}]}]}]], "Input",
 InitializationCell->True,
 CellChangeTimes->{{3.780508949010338*^9, 
  3.7805089945681257`*^9}},ExpressionUUID->"61fb111b-de3e-4add-b4b8-\
87dffb397d30"],

Cell[BoxData[
 RowBox[{"The", " ", "modal", " ", "age", " ", "curve", " ", "is", " ", 
  "plotted", " ", "together", " ", "with", " ", "the", " ", 
  FractionBox[
   RowBox[{"\[Theta]", "-", "\[Theta]1"}], "A"], " ", "curve"}]], "Text",
 CellChangeTimes->{{3.780562271371315*^9, 
  3.780562307784129*^9}},ExpressionUUID->"e363b708-a324-48bb-8d10-\
29f2d1780265"],

Cell[CellGroupData[{

Cell[BoxData[
 RowBox[{"Manipulate", "[", 
  RowBox[{
   RowBox[{"Plot", "[", 
    RowBox[{
     RowBox[{"{", 
      RowBox[{
       RowBox[{
        RowBox[{"ModalAge", "[", "\[Theta]", "]"}], " ", "//.", 
        "\[VeryThinSpace]", 
        RowBox[{"{", "\[VeryThinSpace]", 
         RowBox[{
          RowBox[{"\[Theta]1", "\[Rule]", 
           RowBox[{"\[Theta]0", "+", " ", 
            FractionBox["A", "\[CapitalLambda]"]}]}], ",", 
          RowBox[{"\[Alpha]", "\[Rule]", 
           FractionBox["1", "13"]}], ",", 
          RowBox[{"\[Beta]", "\[Rule]", "b"}], ",", 
          RowBox[{"A", "\[Rule]", "a"}], ",", 
          RowBox[{"\[CapitalLambda]", "\[Rule]", "\[Lambda]"}], ",", 
          RowBox[{"\[Theta]0", "\[Rule]", "\[CurlyRho]"}]}], "}"}]}], ",", 
       RowBox[{
        FractionBox[
         RowBox[{"\[Theta]", "-", "\[Theta]1"}], "A"], "//.", 
        "\[VeryThinSpace]", 
        RowBox[{"{", "\[VeryThinSpace]", 
         RowBox[{
          RowBox[{"\[Theta]1", "\[Rule]", 
           RowBox[{"\[CurlyRho]", "+", " ", 
            FractionBox["A", "\[CapitalLambda]"]}]}], ",", 
          RowBox[{"A", "\[Rule]", "a"}], ",", 
          RowBox[{"\[CapitalLambda]", "\[Rule]", "\[Lambda]"}]}], "}"}]}]}], 
      "}"}], ",", 
     RowBox[{"{", 
      RowBox[{"\[Theta]", ",", "\[CurlyRho]", ",", "440"}], "}"}], ",", 
     RowBox[{"PlotRange", "\[Rule]", "All"}]}], "]"}], ",", 
   RowBox[{"{", 
    RowBox[{
     RowBox[{"{", 
      RowBox[{"a", ",", "1", ",", "\"\<A\>\""}], "}"}], ",", "0.5", ",", 
     "1.5", ",", 
     RowBox[{"Appearance", "\[Rule]", "\"\<Labeled\>\""}]}], "}"}], ",", 
   RowBox[{"{", 
    RowBox[{
     RowBox[{"{", 
      RowBox[{"b", ",", "0.365", ",", "\"\<\[Beta]\>\""}], "}"}], ",", "0.1", 
     ",", "0.6", ",", 
     RowBox[{"Appearance", "\[Rule]", "\"\<Labeled\>\""}]}], "}"}], ",", 
   RowBox[{"{", 
    RowBox[{
     RowBox[{"{", 
      RowBox[{"\[Lambda]", ",", "0.1", ",", "\"\<\[CapitalLambda]\>\""}], 
      "}"}], ",", "0.1", ",", "0.2", ",", 
     RowBox[{"Appearance", "\[Rule]", "\"\<Labeled\>\""}]}], "}"}], ",", 
   RowBox[{"{", 
    RowBox[{
     RowBox[{"{", 
      RowBox[{"\[CurlyRho]", ",", "360", ",", "\"\<\[Theta]0\>\""}], "}"}], 
     ",", "350", ",", "370", ",", 
     RowBox[{"Appearance", "\[Rule]", "\"\<Labeled\>\""}]}], "}"}], ",", 
   RowBox[{"SaveDefinitions", "\[Rule]", "True"}]}], "]"}]], "Input",
 InitializationCell->True,
 CellChangeTimes->{{3.7767075459138055`*^9, 3.776707675218894*^9}, {
   3.7767077065539846`*^9, 3.776707745141926*^9}, {3.7767984387513466`*^9, 
   3.7767985013495703`*^9}, {3.776798914981083*^9, 3.7767989397868867`*^9}, {
   3.7767994793604794`*^9, 3.7767994846024637`*^9}, {3.77680349162313*^9, 
   3.7768034955558667`*^9}, {3.780508246136828*^9, 3.7805082799636307`*^9}, {
   3.7805083255911484`*^9, 3.780508346622979*^9}, 3.7805091048079233`*^9, {
   3.7805598021726923`*^9, 
   3.780559816877124*^9}},ExpressionUUID->"5421ff9c-a63b-4ed2-bdcc-\
39716f2c5f21"],

Cell[BoxData[
 TagBox[
  StyleBox[
   DynamicModuleBox[{$CellContext`a$$ = 0.912, $CellContext`b$$ = 
    0.544, $CellContext`\[Lambda]$$ = 0.1264, $CellContext`\[CurlyRho]$$ = 
    360.65, Typeset`show$$ = True, Typeset`bookmarkList$$ = {}, 
    Typeset`bookmarkMode$$ = "Menu", Typeset`animator$$, Typeset`animvar$$ = 
    1, Typeset`name$$ = "\"untitled\"", Typeset`specs$$ = {{{
       Hold[$CellContext`a$$], 1, "A"}, 0.5, 1.5}, {{
       Hold[$CellContext`b$$], 0.365, "\[Beta]"}, 0.1, 0.6}, {{
       Hold[$CellContext`\[Lambda]$$], 0.1, "\[CapitalLambda]"}, 0.1, 0.2}, {{
       Hold[$CellContext`\[CurlyRho]$$], 360, "\[Theta]0"}, 350, 370}}, 
    Typeset`size$$ = {360., {106., 112.}}, Typeset`update$$ = 0, 
    Typeset`initDone$$, Typeset`skipInitDone$$ = 
    False, $CellContext`a$19265$$ = 0, $CellContext`b$19266$$ = 
    0, $CellContext`\[Lambda]$19267$$ = 0, $CellContext`\[CurlyRho]$19268$$ = 
    0}, 
    DynamicBox[Manipulate`ManipulateBoxes[
     1, StandardForm, 
      "Variables" :> {$CellContext`a$$ = 1, $CellContext`b$$ = 
        0.365, $CellContext`\[Lambda]$$ = 0.1, $CellContext`\[CurlyRho]$$ = 
        360}, "ControllerVariables" :> {
        Hold[$CellContext`a$$, $CellContext`a$19265$$, 0], 
        Hold[$CellContext`b$$, $CellContext`b$19266$$, 0], 
        Hold[$CellContext`\[Lambda]$$, $CellContext`\[Lambda]$19267$$, 0], 
        Hold[$CellContext`\[CurlyRho]$$, $CellContext`\[CurlyRho]$19268$$, 
         0]}, "OtherVariables" :> {
       Typeset`show$$, Typeset`bookmarkList$$, Typeset`bookmarkMode$$, 
        Typeset`animator$$, Typeset`animvar$$, Typeset`name$$, 
        Typeset`specs$$, Typeset`size$$, Typeset`update$$, Typeset`initDone$$,
         Typeset`skipInitDone$$}, "Body" :> Plot[{
         ReplaceRepeated[
          $CellContext`ModalAge[$CellContext`\[Theta]], \
{$CellContext`\[Theta]1 -> $CellContext`\[Theta]0 + \
$CellContext`A/$CellContext`\[CapitalLambda], $CellContext`\[Alpha] -> 
           1/13, $CellContext`\[Beta] -> $CellContext`b$$, $CellContext`A -> \
$CellContext`a$$, $CellContext`\[CapitalLambda] -> $CellContext`\[Lambda]$$, \
$CellContext`\[Theta]0 -> $CellContext`\[CurlyRho]$$}], 
         ReplaceRepeated[($CellContext`\[Theta] - \
$CellContext`\[Theta]1)/$CellContext`A, {$CellContext`\[Theta]1 -> \
$CellContext`\[CurlyRho]$$ + $CellContext`A/$CellContext`\[CapitalLambda], \
$CellContext`A -> $CellContext`a$$, $CellContext`\[CapitalLambda] -> \
$CellContext`\[Lambda]$$}]}, {$CellContext`\[Theta], \
$CellContext`\[CurlyRho]$$, 440}, PlotRange -> All], 
      "Specifications" :> {{{$CellContext`a$$, 1, "A"}, 0.5, 1.5, Appearance -> 
         "Labeled"}, {{$CellContext`b$$, 0.365, "\[Beta]"}, 0.1, 0.6, 
         Appearance -> 
         "Labeled"}, {{$CellContext`\[Lambda]$$, 0.1, "\[CapitalLambda]"}, 
         0.1, 0.2, Appearance -> 
         "Labeled"}, {{$CellContext`\[CurlyRho]$$, 360, "\[Theta]0"}, 350, 
         370, Appearance -> "Labeled"}}, "Options" :> {}, 
      "DefaultOptions" :> {}],
     ImageSizeCache->{411., {195., 201.}},
     SingleEvaluation->True],
    Deinitialization:>None,
    DynamicModuleValues:>{},
    Initialization:>({$CellContext`ModalAge[
         Pattern[$CellContext`\[Theta], 
          Blank[]]] := 
       Max[($CellContext`\[Theta] - $CellContext`\[Theta]1)/$CellContext`A, 
          0] + Max[(1/$CellContext`\[CapitalLambda]) 
          Log[($CellContext`\[Beta] 
             Min[$CellContext`A, $CellContext`\[CapitalLambda] ($CellContext`\
\[Theta] - $CellContext`\[Theta]0)])/($CellContext`\[Alpha] + $CellContext`\
\[CapitalLambda])], 0]}; Typeset`initDone$$ = True),
    SynchronousInitialization->True,
    UndoTrackedVariables:>{Typeset`show$$, Typeset`bookmarkMode$$},
    UnsavedVariables:>{Typeset`initDone$$},
    UntrackedVariables:>{Typeset`size$$}], "Manipulate",
   Deployed->True,
   StripOnInput->False],
  Manipulate`InterpretManipulate[1]]], "Output",
 CellChangeTimes->{
  3.7768035055251646`*^9, 3.780317161625598*^9, 3.7804911542911577`*^9, 
   3.780508005502987*^9, 3.7805082813535233`*^9, 3.7805083489665036`*^9, 
   3.7805085947678123`*^9, 3.7805088169979215`*^9, 3.7805090034596834`*^9, 
   3.780509107495241*^9, {3.7805598035400467`*^9, 3.7805598185488997`*^9}, 
   3.7805607627206583`*^9},ExpressionUUID->"65c5814d-8b85-4446-9176-\
194c23e219c9"]
}, Open  ]]
}, Open  ]]
}, Open  ]]
},
WindowSize->{1127, 689},
WindowMargins->{{131, Automatic}, {23, Automatic}},
FrontEndVersion->"11.2 for Microsoft Windows (64-bit) (September 10, 2017)",
StyleDefinitions->FrontEnd`FileName[{"Report"}, "StandardReport.nb", 
  CharacterEncoding -> "UTF-8"]
]
(* End of Notebook Content *)

(* Internal cache information *)
(*CellTagsOutline
CellTagsIndex->{}
*)
(*CellTagsIndex
CellTagsIndex->{}
*)
(*NotebookFileOutline
Notebook[{
Cell[CellGroupData[{
Cell[1510, 35, 299, 7, 155, "Title",ExpressionUUID->"4b842629-12ca-4fb3-9916-f4cd004afc45"],
Cell[1812, 44, 247, 6, 65, "Subtitle",ExpressionUUID->"208b26fc-afe1-46d7-a92d-abe32610fef1"],
Cell[2062, 52, 285, 6, 33, "Text",ExpressionUUID->"c9cf37f4-8fc4-41b6-bb80-29a53e45be7a"],
Cell[2350, 60, 582, 10, 77, "Text",ExpressionUUID->"0064d11b-757f-4896-93f6-db41f8147880"],
Cell[CellGroupData[{
Cell[2957, 74, 185, 3, 67, "Chapter",ExpressionUUID->"00369a74-86ec-41e3-88b7-4a6f06ad123d"],
Cell[3145, 79, 3594, 78, 155, "Text",ExpressionUUID->"7781d49e-ea5e-4269-b06f-d6bef7d370d1"],
Cell[6742, 159, 2527, 67, 95, "Input",ExpressionUUID->"3eeefdd3-c614-4e2f-a9bc-ffaff7edafbb",
 InitializationCell->True]
}, Open  ]],
Cell[CellGroupData[{
Cell[9306, 231, 157, 3, 67, "Chapter",ExpressionUUID->"80f2d523-6d0b-4d7b-a5d6-971d478aafa5"],
Cell[CellGroupData[{
Cell[9488, 238, 241, 4, 69, "Section",ExpressionUUID->"350606bb-2e08-499d-94fa-d090cfbd0cb0"],
Cell[CellGroupData[{
Cell[9754, 246, 3390, 73, 132, "Input",ExpressionUUID->"d5610d9c-1d00-4816-8a21-020cc0bc0608",
 InitializationCell->True],
Cell[13147, 321, 5159, 92, 446, "Output",ExpressionUUID->"70e2c71a-58e6-448b-ba15-5b0190dd85ea"]
}, Open  ]]
}, Open  ]],
Cell[CellGroupData[{
Cell[18355, 419, 184, 3, 69, "Section",ExpressionUUID->"6b94b792-a36b-4cd3-bc80-6ad9d9535a24"],
Cell[CellGroupData[{
Cell[18564, 426, 3449, 73, 132, "Input",ExpressionUUID->"275c7960-e4da-44b9-ba32-b10e12500197",
 InitializationCell->True],
Cell[22016, 501, 4980, 91, 446, "Output",ExpressionUUID->"2c66fc81-5a79-4299-9919-96282d97e0f6"]
}, Open  ]],
Cell[27011, 595, 6834, 162, 284, "Text",ExpressionUUID->"464b4d9a-53aa-4de3-bc4f-997cb4be354b"]
}, Open  ]]
}, Open  ]],
Cell[CellGroupData[{
Cell[33894, 763, 239, 4, 67, "Chapter",ExpressionUUID->"b09a3313-ec68-4def-8a0b-43b2c87b1621"],
Cell[CellGroupData[{
Cell[34158, 771, 3143, 68, 113, "Input",ExpressionUUID->"7e669c3b-775b-411b-b444-cc600741cad6",
 InitializationCell->True],
Cell[37304, 841, 4771, 88, 560, "Output",ExpressionUUID->"3cdf87bb-0244-4bb2-a9d9-6d02ad90c767"]
}, Open  ]]
}, Open  ]],
Cell[CellGroupData[{
Cell[42124, 935, 406, 6, 67, "Chapter",ExpressionUUID->"65889081-75c0-4f05-ae9d-3ced8225cccc"],
Cell[CellGroupData[{
Cell[42555, 945, 246, 4, 69, "Section",ExpressionUUID->"49a9ea81-f79e-4cbc-b27a-1cc93ffdabd7"],
Cell[CellGroupData[{
Cell[42826, 953, 831, 22, 39, "Input",ExpressionUUID->"335a8e78-088d-4042-af62-4d1b067c96b7"],
Cell[43660, 977, 932, 29, 81, "Output",ExpressionUUID->"0939f8c4-a8a7-48cd-b9d9-6e8d83ee3ba2"]
}, Open  ]]
}, Open  ]],
Cell[CellGroupData[{
Cell[44641, 1012, 268, 4, 69, "Section",ExpressionUUID->"23015b51-d638-4ed6-bfbd-4f92d0bb9316"],
Cell[CellGroupData[{
Cell[44934, 1020, 392, 6, 33, "Subsection",ExpressionUUID->"c6ee4b5d-65c5-455b-8e78-dfd3ac080cfc"],
Cell[45329, 1028, 198, 3, 33, "Text",ExpressionUUID->"fdd5383f-495b-49e6-aa94-6c9e33ed4faa"],
Cell[CellGroupData[{
Cell[45552, 1035, 1603, 43, 62, "Input",ExpressionUUID->"45788ab6-b913-4357-bb5e-93b46866de11"],
Cell[47158, 1080, 1230, 35, 63, "Output",ExpressionUUID->"0a2ad285-f54e-4da7-afb1-c74f89f642d4"]
}, Open  ]],
Cell[48403, 1118, 683, 16, 39, "Text",ExpressionUUID->"50496439-cb9f-48c2-8f54-325085d4afa0"]
}, Open  ]],
Cell[CellGroupData[{
Cell[49123, 1139, 421, 6, 33, "Subsection",ExpressionUUID->"4f709d11-8991-4045-8acf-83ddc3a763f6"],
Cell[CellGroupData[{
Cell[49569, 1149, 1114, 27, 62, "Input",ExpressionUUID->"fedcd1fb-8b64-467e-8ec0-001ed6ab0fa0"],
Cell[50686, 1178, 780, 25, 90, "Output",ExpressionUUID->"7d22e2b1-8524-45e3-8424-0ea61d1181a5"]
}, Open  ]],
Cell[51481, 1206, 1197, 26, 68, "Text",ExpressionUUID->"784746ff-374c-44af-94b6-8b0ad1f03400"]
}, Open  ]]
}, Open  ]],
Cell[CellGroupData[{
Cell[52727, 1238, 342, 5, 69, "Section",ExpressionUUID->"869de240-e246-4f88-8a84-8b69843583c3"],
Cell[CellGroupData[{
Cell[53094, 1247, 836, 23, 39, "Input",ExpressionUUID->"6a8404d0-4df9-4fda-8b62-338f5ac43052"],
Cell[53933, 1272, 441, 12, 39, "Output",ExpressionUUID->"5fa078d3-7cfd-4d73-b13e-38dae4e7acc2"]
}, Open  ]],
Cell[CellGroupData[{
Cell[54411, 1289, 223, 5, 33, "Subsection",ExpressionUUID->"8ddc7c60-240a-4a87-8251-9a209e283a4a"],
Cell[54637, 1296, 2611, 76, 210, "Input",ExpressionUUID->"c7696fcc-cdb5-4ba2-b455-2b9760e13587",
 InitializationCell->True],
Cell[57251, 1374, 1962, 53, 93, "Text",ExpressionUUID->"86ce030f-ef92-4766-84a1-fc817d163ccf"]
}, Open  ]]
}, Open  ]]
}, Open  ]],
Cell[CellGroupData[{
Cell[59274, 1434, 221, 3, 67, "Chapter",ExpressionUUID->"7fad5e09-8f8e-474b-9c45-f7396711b77f"],
Cell[CellGroupData[{
Cell[59520, 1441, 172, 3, 69, "Section",ExpressionUUID->"5efd3c13-e37e-4e02-9547-65f44ca71744"],
Cell[59695, 1446, 4388, 121, 186, "Input",ExpressionUUID->"62a3ecfb-ddaa-40c4-bfaf-e3a256be4af9",
 InitializationCell->True]
}, Open  ]],
Cell[CellGroupData[{
Cell[64120, 1572, 163, 3, 69, "Section",ExpressionUUID->"2c660588-a942-43e3-a1da-1b64c15827aa"],
Cell[CellGroupData[{
Cell[64308, 1579, 3179, 69, 113, "Input",ExpressionUUID->"554bcbeb-5eca-435c-be63-b02d5784b0f6",
 InitializationCell->True],
Cell[67490, 1650, 3872, 70, 560, "Output",ExpressionUUID->"c1656431-dca3-493f-a989-bd33d7f5197a"]
}, Open  ]]
}, Open  ]],
Cell[CellGroupData[{
Cell[71411, 1726, 184, 3, 69, "Section",ExpressionUUID->"2cc81282-9da0-4885-9c78-9cf60985c362"],
Cell[CellGroupData[{
Cell[71620, 1733, 3141, 69, 132, "Input",ExpressionUUID->"6440052f-01b3-4211-94b9-5df0a3a02e3b",
 InitializationCell->True],
Cell[74764, 1804, 6301, 114, 450, "Output",ExpressionUUID->"4807522d-91a0-49b9-8501-843b5ba2b458"]
}, Open  ]]
}, Open  ]]
}, Open  ]],
Cell[CellGroupData[{
Cell[81126, 1925, 404, 6, 67, "Chapter",ExpressionUUID->"a78d3a04-98f5-4f39-ae53-ce7070c88d9f"],
Cell[CellGroupData[{
Cell[81555, 1935, 246, 4, 69, "Section",ExpressionUUID->"811db1ff-74a3-4036-9db6-315c487ccdfb"],
Cell[CellGroupData[{
Cell[81826, 1943, 824, 25, 39, "Input",ExpressionUUID->"6e1dcdd2-0610-4537-a3e2-a00b68e25821"],
Cell[82653, 1970, 898, 26, 79, "Output",ExpressionUUID->"737742ed-bc19-41f5-8f82-1fa39ef9f81b"]
}, Open  ]]
}, Open  ]],
Cell[CellGroupData[{
Cell[83600, 2002, 268, 4, 69, "Section",ExpressionUUID->"328d4a86-1a09-44c6-a226-922c8750058b"],
Cell[CellGroupData[{
Cell[83893, 2010, 333, 5, 33, "Subsection",ExpressionUUID->"200a26cf-4b03-4423-b796-5092098a7d9b"],
Cell[CellGroupData[{
Cell[84251, 2019, 1623, 44, 62, "Input",ExpressionUUID->"ff942149-91c0-490f-9fb0-9bb9cf49422d"],
Cell[85877, 2065, 1211, 32, 86, "Output",ExpressionUUID->"0442df89-9f9a-4ece-8c6a-24c9ec58d648"]
}, Open  ]]
}, Open  ]],
Cell[CellGroupData[{
Cell[87137, 2103, 364, 5, 33, "Subsection",ExpressionUUID->"2ffbfaa1-afd7-438c-be37-6ce62fad66f4"],
Cell[CellGroupData[{
Cell[87526, 2112, 1661, 48, 62, "Input",ExpressionUUID->"e5d69a03-dec6-4fb2-b9f0-98abeddef99e"],
Cell[89190, 2162, 1386, 40, 106, "Output",ExpressionUUID->"7a5a5d94-1f3f-47ac-9c5b-4b7f16cb494e"]
}, Open  ]],
Cell[90591, 2205, 1260, 36, 95, "Text",ExpressionUUID->"c000d6a2-7ed1-4c4f-9beb-b6b9af979965"]
}, Open  ]],
Cell[CellGroupData[{
Cell[91888, 2246, 209, 4, 33, "Subsection",ExpressionUUID->"f06ec256-c77d-4014-927e-58faa4334509"],
Cell[CellGroupData[{
Cell[92122, 2254, 205, 3, 32, "Subsubsection",ExpressionUUID->"57b2f3e2-f8c1-4c57-bb83-83e683d5c1ca"],
Cell[92330, 2259, 3353, 99, 269, "Input",ExpressionUUID->"9d0b4dfa-f0ea-4362-ace5-51e3fe15721d",
 InitializationCell->True]
}, Closed]],
Cell[CellGroupData[{
Cell[95720, 2363, 181, 3, 32, "Subsubsection",ExpressionUUID->"754de79c-2935-4a6c-bb47-3f088351eb09"],
Cell[95904, 2368, 1102, 22, 76, "Input",ExpressionUUID->"66d5669c-5161-4653-919f-eee583243123"]
}, Open  ]],
Cell[CellGroupData[{
Cell[97043, 2395, 187, 3, 32, "Subsubsection",ExpressionUUID->"9286bf45-5db9-4cd0-91dd-8c6e1c3f06bd"],
Cell[97233, 2400, 709, 15, 74, "Text",ExpressionUUID->"45b6e362-0c36-463b-959a-bf15cb35087a"],
Cell[CellGroupData[{
Cell[97967, 2419, 8383, 239, 678, "Input",ExpressionUUID->"25b6c3b0-f5a8-423d-99a8-cd5c7159ddb9",
 InitializationCell->True],
Cell[106353, 2660, 6989, 121, 422, "Output",ExpressionUUID->"010f8b2f-f465-4c3c-9110-56bdb7d4d6b9"]
}, Open  ]],
Cell[113357, 2784, 5891, 123, 154, "Text",ExpressionUUID->"ad386b5b-589f-4589-9096-21a4a2415beb"]
}, Open  ]],
Cell[CellGroupData[{
Cell[119285, 2912, 198, 3, 32, "Subsubsection",ExpressionUUID->"fd34bb06-f856-4e0a-ad7c-7c5ab2a0670c"],
Cell[119486, 2917, 195, 3, 27, "Item",ExpressionUUID->"ff06dca1-14fd-48b3-9b7c-c341435f7aae"],
Cell[CellGroupData[{
Cell[119706, 2924, 3531, 79, 156, "Input",ExpressionUUID->"01d9c5e4-dfb2-4a5e-928c-8b4739a5c827",
 InitializationCell->True],
Cell[123240, 3005, 4675, 84, 70, "Output",ExpressionUUID->"98033c26-1d37-4dc2-bc70-3c45609ac9fe"]
}, Open  ]],
Cell[127930, 3092, 284, 7, 33, "Text",ExpressionUUID->"1aa62362-fd33-4017-b329-385bb23196c0"],
Cell[128217, 3101, 168, 3, 27, "Item",ExpressionUUID->"d7377737-0593-4406-a3eb-301c07e09aeb"],
Cell[CellGroupData[{
Cell[128410, 3108, 4322, 95, 158, "Input",ExpressionUUID->"9b2eee7b-58d7-4300-804b-0553335bd6d6",
 InitializationCell->True],
Cell[132735, 3205, 5847, 101, 70, "Output",ExpressionUUID->"92a14862-0dd7-4d3b-bbf5-2b1285e6dfb5"]
}, Open  ]],
Cell[138597, 3309, 270, 6, 33, "Text",ExpressionUUID->"5758f0ed-1d54-4b60-831b-c4eb491fd5ff"]
}, Closed]]
}, Open  ]]
}, Open  ]]
}, Open  ]],
Cell[CellGroupData[{
Cell[138940, 3323, 459, 7, 67, "Chapter",ExpressionUUID->"b65f26ef-55e2-4a84-ac82-b23e4dc95dde"],
Cell[139402, 3332, 476, 8, 33, "Text",ExpressionUUID->"b06bcadf-1c7e-49f6-b81a-6b0e8929ebb1"],
Cell[139881, 3342, 858, 25, 63, "Input",ExpressionUUID->"61fb111b-de3e-4add-b4b8-87dffb397d30",
 InitializationCell->True],
Cell[140742, 3369, 362, 7, 50, "Text",ExpressionUUID->"e363b708-a324-48bb-8d10-29f2d1780265"],
Cell[CellGroupData[{
Cell[141129, 3380, 2983, 70, 123, "Input",ExpressionUUID->"5421ff9c-a63b-4ed2-bdcc-39716f2c5f21",
 InitializationCell->True],
Cell[144115, 3452, 4318, 78, 422, "Output",ExpressionUUID->"65c5814d-8b85-4446-9176-194c23e219c9"]
}, Open  ]]
}, Open  ]]
}, Open  ]]
}
]
*)

(* NotebookSignature wx0TWlAbfB0M7DKv15rRkjMb *)
