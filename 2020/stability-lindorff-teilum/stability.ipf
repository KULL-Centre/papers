#pragma rtGlobals=3		// Use modern global access method and strict wave access.
#include <KBColorizeTraces>
Function H(w,x) : FitFunc
	Wave w
	Variable x
	//CurveFitDialog/ These comments were created by the Curve Fitting dialog. Altering them will
	//CurveFitDialog/ make the function less convenient to work with in the Curve Fitting dialog.
	//CurveFitDialog/ Equation:
	//CurveFitDialog/ f(x) = H_Tm+Cp*(x-Tm)
	//CurveFitDialog/ End of Equation
	//CurveFitDialog/ Independent Variables 1
	//CurveFitDialog/ x
	//CurveFitDialog/ Coefficients 3
	//CurveFitDialog/ w[0] = H_Tm
	//CurveFitDialog/ w[1] = Tm
	//CurveFitDialog/ w[2] = Cp
	return w[0]+w[2]*(x-w[1])
End
Function TS(w,x) : FitFunc
	Wave w
	Variable x
	//CurveFitDialog/ These comments were created by the Curve Fitting dialog. Altering them will
	//CurveFitDialog/ make the function less convenient to work with in the Curve Fitting dialog.
	//CurveFitDialog/ Equation:
	//CurveFitDialog/ f(x) = (H_Tm/Tm+Cp*ln(x/Tm))*x
	//CurveFitDialog/ End of Equation
	//CurveFitDialog/ Independent Variables 1
	//CurveFitDialog/ x
	//CurveFitDialog/ Coefficients 3
	//CurveFitDialog/ w[0] = H_Tm
	//CurveFitDialog/ w[1] = Tm
	//CurveFitDialog/ w[2] = Cp
	return (w[0]/w[1]+w[2]*ln(x/w[1]))*x
End
Function G(w,x,y)
	Wave w
	Variable x
	Variable y
	return H({w[0],w[1],w[2]},x)-TS({w[0],w[1],w[2]},x)-w[3]*y
End
Function f(w,x,y)
	Wave w
	Variable x
	Variable y
	return ((10-0.01*(x-273)-0.05*y)+(3-0.05*(x-273)+0.00025*(x-273)^2-0.1*y)*exp(-G(w,x,y)/8.314e-3/x))/(1+exp(-G(w,x,y)/8.314e-3/x))
End
Window Graph1() : Graph
	PauseUpdate; Silent 1		// building window...
	Display /W=(156,45,442,234) dH_T,dS_T,dG_T
	ModifyGraph gFont="Helvetica",gfSize=10,width={Aspect,1.41},height=140.031
	ModifyGraph lSize=2
	ModifyGraph rgb(dH_T)=(0,0,65535),rgb(dG_T)=(3,52428,1)
	ModifyGraph zero(left)=1
	ModifyGraph mirror=2
	ModifyGraph lblMargin(left)=9
	Label left "\\F'Symbol'Δ\\F'Helvetica'\\Br\\M\\f02G\\f00, \\F'Symbol'Δ\\F'Helvetica'\\Br\\M\\f02H\\f00 or T\\F'Symbol'Δ\\F'Helvetica'\\Br\\M\\f02S\\f00 (kJ mol\\S-1\\M)"
	Label bottom "Temperature (K)"
	SetAxis left*,628.72796
	SetAxis bottom*,353
	TextBox/C/N=text0/F=0/A=MC/X=45.18/Y=27.14 "\\K(0,0,65535)\\F'Symbol'Δ\\F'Helvetica'\\Br\\M\\f02H"
	TextBox/C/N=text1/F=0/A=MC/X=36.55/Y=42.14 "\\K(65535,0,0)\\F'Helvetica'T\\F'Symbol'Δ\\F'Helvetica'\\Br\\M\\f02S"
	TextBox/C/N=text2/F=0/A=MC/X=39.59/Y=-18.57 "\\K(3,52427,1)\\F'Symbol'Δ\\F'Helvetica'\\Br\\M\\f02G"
EndMacro
Window Graph1_1() : Graph
	PauseUpdate; Silent 1		// building window...
	Display /W=(467,45,747,234) dH_T,dS_T,dG_T
	ModifyGraph gFont="Helvetica",gfSize=10,width={Aspect,1.41},height=140.031
	ModifyGraph lSize=2
	ModifyGraph rgb(dH_T)=(0,0,65535),rgb(dG_T)=(3,52428,1)
	ModifyGraph zero(left)=1
	ModifyGraph mirror=2
	ModifyGraph lblMargin(left)=9
	Label left "\\F'Symbol'Δ\\F'Helvetica'\\Br\\M\\f02G\\f00, \\F'Symbol'Δ\\F'Helvetica'\\Br\\M\\f02H\\f00 or T\\F'Symbol'Δ\\F'Helvetica'\\Br\\M\\f02S\\f00 (kJ mol\\S-1\\M)"
	Label bottom "Temperature (K)"
	SetAxis left -20,20
	SetAxis bottom*,353
	TextBox/C/N=text0/F=0/A=MC/X=-13.71/Y=-25.00 "\\K(0,0,65535)\\F'Symbol'Δ\\F'Helvetica'\\Br\\M\\f02H"
	TextBox/C/N=text1/F=0/A=MC/X=3.05/Y=-17.86 "\\K(65535,0,0)\\F'Helvetica'T\\F'Symbol'Δ\\F'Helvetica'\\Br\\M\\f02S"
	TextBox/C/N=text2/F=0/A=MC/X=41.62/Y=-17.86 "\\K(3,52427,1)\\F'Symbol'Δ\\F'Helvetica'\\Br\\M\\f02G"
EndMacro
Window Graph5() : Graph
	PauseUpdate; Silent 1		// building window...
	Display /W=(35,45,304,234) f_3d[*][0],f_3d[*][11],f_3d[*][23],f_3d[*][35],f_3d[*][47]
	AppendToGraph f_3d[*][59],f_3d[*][71],f_3d[*][93],f_3d[*][115],f_3d[*][127]
	ModifyGraph gFont="Helvetica",gfSize=10,width={Aspect,1.41},height=140.031
	ModifyGraph lSize=2
	ModifyGraph rgb(f_3d)=(0,60652,60652),rgb(f_3d#1)=(0,41120,63222),rgb(f_3d#2)=(0,65535,0)
	ModifyGraph rgb(f_3d#3)=(0,51400,0),rgb(f_3d#4)=(65535,65535,0),rgb(f_3d#5)=(59367,49344,0)
	ModifyGraph rgb(f_3d#7)=(54998,0,0),rgb(f_3d#8)=(65535,0,65535),rgb(f_3d#9)=(39321,21845,51657)
	ModifyGraph mirror=2
	ModifyGraph lblMargin(left)=9
	Label left "α (arbitrary unit)"
	Label bottom "Temperature (K)"
	SetAxis left 0,10
	SetAxis bottom 273,373
	ColorScale/C/N=text0/F=0/B=1/A=MC/X=32.00/Y=14.00  ctab={0,4,dBZ14,0}
	AppendText "[denaturant] (M)"
EndMacro
Window Graph6() : Graph
	PauseUpdate; Silent 1		// building window...
	Display /W=(311,45,580,234) f_3d[102][*],f_3d[89][*],f_3d[76][*],f_3d[64][*],f_3d[51][*]
	AppendToGraph f_3d[38][*],f_3d[26][*],f_3d[13][*],f_3d[0][*]
	ModifyGraph gFont="Helvetica",gfSize=10,width={Aspect,1.41},height=140.031
	ModifyGraph lSize=2
	ModifyGraph rgb(f_3d)=(40230,7902,14081),rgb(f_3d#1)=(52680,14908,18146),rgb(f_3d#2)=(62127,41599,31195)
	ModifyGraph rgb(f_3d#3)=(63377,53677,39069),rgb(f_3d#4)=(62103,60354,47563),rgb(f_3d#5)=(53044,59346,61212)
	ModifyGraph rgb(f_3d#6)=(15060,47929,60543),rgb(f_3d#7)=(0,39474,56165),rgb(f_3d#8)=(0,24098,44385)
	ModifyGraph mirror=2
	ModifyGraph lblMargin(left)=9
	Label left "α (arbitrary unit)"
	Label bottom "[denaturant] (M)"
	SetAxis left 0,10
	ColorScale/C/N=text0/F=0/B=1/A=MC/X=32.00/Y=15.00
	ColorScale/C/N=text0  ctab={273,353,EOSSpectral11,1}
	AppendText "Temperature (K)"
EndMacro
