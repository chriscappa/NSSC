#pragma rtGlobals=3		// Use modern global access method and strict wave access.
#include <KBColorizeTraces>

// Functions to load and do basic processing on data from a FIGAERO-CIMS
// Specific functions are set up to deal with data provided from Joel Thornton's group (current or former members)
// Initially written by Ziyue Li (2018-19)

// Version updates
// v1.0 - initial functions written

//Run this function only after all the data are loaded in different folders
//!!! Set current folder to the parent folder that contains wave "Names_origin"
//!!! Current folder must have Subfolders entitled "Normal" and "TD"
//This function creates waves of compound information in the folders "Mass" and "Signal" that are in the subfolder "Normal" and "TD"
Function StandardProcess_Basic()
	If(Datafolderexists(":Normal")==0 || Datafolderexists(":TD")==0)	
		print "Please make sure that you are in the right data folder"
		abort
	endif
	If(Datafolderexists(":Normal:Mass")==0 || Datafolderexists(":Normal:Signal")==0 ||Datafolderexists(":TD:Mass")==0 || Datafolderexists(":TD:Signal")==0)	
		print "Please make sure that you have the right sub data folders"
		abort
	endif
	If(Waveexists($("Names_origin"))==0)
		print "Please make sure that you have the right waves Names_origin in the current folder"
		abort
	endif
	
	string parentfolder_Str=Getdatafolder(1)
	
	//Translate names to a matrix containing numbers of elements in every compound/ion
	wave/t Names_origin
	MakeEleNumMatrix_s(Names_origin)
	wave Ele_num
	GetCmpMass_s(Ele_num)  //Calculate molar mass of each compound by removing Iodine from the ion
	GetIonMass_s(Ele_num)  //Calculate molar mass of each ion	
	wave MolarMass_cmp
	wave MolarMass_ion
	
	//Copy the waves in parent folder to subfolders
	string subfolder_str1=":Normal:Mass"
	string subfolder_str2=":Normal:Signal"
	string subfolder_str3=":TD:Mass"
	string subfolder_str4=":TD:Signal"
	string wavestr1=":Names_origin"
	string wavestr2=":Ele_num"
	string wavestr3=":MolarMass_cmp"
	string wavestr4=":MolarMass_ion"
	duplicate/o $(wavestr1) $(subfolder_str1+wavestr1)
	duplicate/o $(wavestr2) $(subfolder_str1+wavestr2)
	duplicate/o $(wavestr3) $(subfolder_str1+wavestr3)
	duplicate/o $(wavestr4) $(subfolder_str1+wavestr4)
	duplicate/o $(wavestr1) $(subfolder_str2+wavestr1)
	duplicate/o $(wavestr2) $(subfolder_str2+wavestr2)
	duplicate/o $(wavestr3) $(subfolder_str2+wavestr3)
	duplicate/o $(wavestr4) $(subfolder_str2+wavestr4)
	duplicate/o $(wavestr1) $(subfolder_str3+wavestr1)
	duplicate/o $(wavestr2) $(subfolder_str3+wavestr2)
	duplicate/o $(wavestr3) $(subfolder_str3+wavestr3)
	duplicate/o $(wavestr4) $(subfolder_str3+wavestr4)
	duplicate/o $(wavestr1) $(subfolder_str4+wavestr1)
	duplicate/o $(wavestr2) $(subfolder_str4+wavestr2)
	duplicate/o $(wavestr3) $(subfolder_str4+wavestr3)
	duplicate/o $(wavestr4) $(subfolder_str4+wavestr4)
	
	//Analyze data in datafolder :Normal:Mass first
	Setdatafolder subfolder_str1
	CalculateHC() //calculate HC_ratio and other ratios
	Setdatafolder parentfolder_str
	Setdatafolder subfolder_str2
	CalculateHC()
	Setdatafolder parentfolder_str
	Setdatafolder subfolder_str3
	CalculateHC()
	Setdatafolder parentfolder_str
	Setdatafolder subfolder_str4
	CalculateHC()
	
	Setdatafolder parentfolder_str
	print "Correct!"
	
End


//Run this function only after all the data are loaded in different folders
//!!! Set current folder to the parent folder that contains wave "Names_origin"
//!!! Current folder must have Subfolders entitled "Normal" and "TD"
//This function analyze data in the folders "Mass" and "Signal" that are in the subfolder "Normal" and "TD"
Function StandardProcess_Old()
	If(Datafolderexists(":Normal")==0 || Datafolderexists(":TD")==0)	
		print "Please make sure that you are in the right data folder"
		abort
	endif
	If(Datafolderexists(":Normal:Mass")==0 || Datafolderexists(":Normal:Signal")==0 ||Datafolderexists(":TD:Mass")==0 || Datafolderexists(":TD:Signal")==0)	
		print "Please make sure that you have the right sub data folders"
		abort
	endif
	If(Waveexists($("Names_origin"))==0)
		print "Please make sure that you have the right waves Names_origin in the current folder"
		abort
	endif
	
	string parentfolder_Str=Getdatafolder(1)
	
	//Translate names to a matrix containing numbers of elements in every compound/ion
	wave/t Names_origin
	MakeEleNumMatrix_s(Names_origin)
	wave Ele_num
	GetCmpMass_s(Ele_num)  //Calculate molar mass of each compound by removing Iodine from the ion
	GetIonMass_s(Ele_num)  //Calculate molar mass of each ion	
	wave MolarMass_cmp
	wave MolarMass_ion
	
	//Copy the waves in parent folder to subfolders
	string subfolder_str1=":Normal:Mass"
	string subfolder_str2=":Normal:Signal"
	string subfolder_str3=":TD:Mass"
	string subfolder_str4=":TD:Signal"
	string wavestr1=":Names_origin"
	string wavestr2=":Ele_num"
	string wavestr3=":MolarMass_cmp"
	string wavestr4=":MolarMass_ion"
	duplicate/o $(wavestr1) $(subfolder_str1+wavestr1)
	duplicate/o $(wavestr2) $(subfolder_str1+wavestr2)
	duplicate/o $(wavestr3) $(subfolder_str1+wavestr3)
	duplicate/o $(wavestr4) $(subfolder_str1+wavestr4)
	duplicate/o $(wavestr1) $(subfolder_str2+wavestr1)
	duplicate/o $(wavestr2) $(subfolder_str2+wavestr2)
	duplicate/o $(wavestr3) $(subfolder_str2+wavestr3)
	duplicate/o $(wavestr4) $(subfolder_str2+wavestr4)
	duplicate/o $(wavestr1) $(subfolder_str3+wavestr1)
	duplicate/o $(wavestr2) $(subfolder_str3+wavestr2)
	duplicate/o $(wavestr3) $(subfolder_str3+wavestr3)
	duplicate/o $(wavestr4) $(subfolder_str3+wavestr4)
	duplicate/o $(wavestr1) $(subfolder_str4+wavestr1)
	duplicate/o $(wavestr2) $(subfolder_str4+wavestr2)
	duplicate/o $(wavestr3) $(subfolder_str4+wavestr3)
	duplicate/o $(wavestr4) $(subfolder_str4+wavestr4)
	
	//Analyze data in datafolder :Normal:Mass first
	Setdatafolder subfolder_str1
	CalculateHC() //calculate HC_ratio and other ratios
	wave OC_ratio
	wave HC_ratio
	wave NC_ratio
	wave OSc
	wave time_s
	wave temp_c
	wave cmp_matrix
	if(waveexists(time_s)==0 || waveexists(temp_c)==0 || waveexists(cmp_matrix)==0)
		print "Please make sure that subfolder "+subfolder_str1+"  has the required waves"
		abort
	endif
	IntegrateEachCol1(cmp_matrix,time_s)  //integrate mass of each ion
	wave cmp_mass
	wave cmp_MFR
	wave cmp_matrix_INT
	display cmp_mass vs molarmass_cmp
	ModifyGraph mirror=1,fSize=14,standoff=0;DelayUpdate
	Label left "Mass, ng";DelayUpdate
	Label bottom "m/z";DelayUpdate
	SetAxis left 0,*
	ModifyGraph mode=1
	wavestats/q cmp_mass
	duplicate/o cmp_mass cmp_mass_norm
	wave cmp_mass_norm
	cmp_mass_norm/=V_sum
	FilterbyFraction(cmp_mass_norm,0.005,molarmass_ion)
	FilterbyFraction(cmp_mass_norm,0.005,molarmass_cmp)
	wave cmp_mass_norm_filter
	wave molarmass_ion_filter
	wave molarmass_cmp_filter
	
	
	Setdatafolder parentfolder_str
	print "Correct!"
	
End


Function MakeNewFolders()
	Setdatafolder root:
	Newdatafolder/o Normal
	Newdatafolder/o TD
	Setdatafolder root:Normal:
	Newdatafolder/o Signal
	Newdatafolder/o Mass
	Setdatafolder root:Normal:Signal:
	CreateNewFolders_signal()
	Setdatafolder root:Normal:Mass:
	CreateNewFolders_mass()
	Setdatafolder root:TD:
	Newdatafolder/o Signal
	Newdatafolder/o Mass
	Setdatafolder root:TD:Signal:
	CreateNewFolders_signal()
	Setdatafolder root:TD:Mass:
	CreateNewFolders_mass()
End

Function CreateNewFolders_signal()
	NewDataFolder/o $("RawSignal")
	NewDataFolder/o $("RawSignal_bgd")
	NewDataFolder/o $("Signal")
End


Function CreateNewFolders_mass()
	NewDataFolder/o $("RawMass")
	NewDataFolder/o $("RawMass_bgd")
	NewDataFolder/o $("Mass")
End


//For Sigi's data
Function CreateNewFolders_s()
	NewDataFolder/o $("RawMass")
	NewDataFolder/o $("RawMass_norm")
	NewDataFolder/o $("RawSignal")
	NewDataFolder/o $("Mass")
	NewDataFolder/o $("Mass_norm")
	NewDataFolder/o $("Signal")
End

// This function finds the number of a specific element in the compound
//c_str is the elemental composition of the compound
//e_str is the targeted element, e.g. C, H, O, I etc. 
//This function returns e_num which is the number of element e_str in the compound c_str
Function ElementNumber(c_str,e_str)
	string c_str
	string e_str
	
	variable c_pos,e_num
	variable len
	len=strlen(c_str)
	c_pos=strsearch(c_str,e_str,0,2)
	if(c_pos==-1)
		e_num=0
	else
		string substring
		substring=c_str[c_pos+1,len-1]
		e_num=str2num(substring)
		if(numtype(e_num)==2)
			e_num=1
		endif
	endif
	return e_num	
End

//Make a matrix that stores the number of each element, C, H, O, N, I of all the compounds
Function MakeEleNumMatrix(Names_origin)
	wave/t Names_origin
	variable wavelen=numpnts(names_origin)
	
	make/o/d/n=(wavelen,5) Ele_num
	string cmp_str
	variable i
	for(i=0;i<wavelen;i+=1)
		cmp_str=names_origin[i]
		Ele_num[i][0]=ElementNumber(cmp_str,"C")
		Ele_num[i][1]=ElementNumber(cmp_str,"H")
		Ele_num[i][2]=ElementNumber(cmp_str,"O")
		Ele_num[i][3]=ElementNumber(cmp_str,"N")
		Ele_num[i][4]=ElementNumber(cmp_str,"I")
	endfor
End

//Make a matrix that stores the number of each element, C, H, O, N, I,S of all the compounds
Function MakeEleNumMatrix_s(Names_origin)
	wave/t Names_origin
	variable wavelen=numpnts(names_origin)
	
	make/o/d/n=(wavelen,6) Ele_num
	string cmp_str
	variable i
	for(i=0;i<wavelen;i+=1)
		cmp_str=names_origin[i]
		Ele_num[i][0]=ElementNumber(cmp_str,"C")
		Ele_num[i][1]=ElementNumber(cmp_str,"H")
		Ele_num[i][2]=ElementNumber(cmp_str,"O")
		Ele_num[i][3]=ElementNumber(cmp_str,"N")
		Ele_num[i][4]=ElementNumber(cmp_str,"I")
		Ele_num[i][5]=ElementNumber(cmp_str,"S")
	endfor
End

Function SigiCmpFilter()
	wave CmpN_all
	wave/t Names_all
	wave ele_num_all
	
	duplicate/o CmpN_all CmpN_origin
	duplicate/o/t Names_all Names_origin
	duplicate/o Ele_num_all Ele_num
	
	variable len=numpnts(cmpn_all)
	variable i, count
	count=0
	string tempstr
	
	for(i=0;i<len;i+=1)
		if(ele_num_all[i][0]!=0 && ele_num_all[i][1]!=0 && ele_num_all[i][2]!=0 && ele_num_all[i][4]==1)
			cmpn_origin[count]=cmpn_all[i]
			tempstr=names_all[i]
			names_origin[count]=tempstr
			ele_num[count][]=ele_num_all[i][q]
			count+=1
		endif
	endfor
	
	Deletepoints count,(len-count),cmpn_origin
	Deletepoints count,(len-count),names_origin
	Deletepoints/m=0 count,(len-count),ele_num
	
End

//Calculate the molar mass of each compound
Function GetIonMass(Ele_num)
	wave Ele_num
	variable num_cmp=dimsize(ele_num,0)
	make/o/d/n=(num_cmp) MolarMass_ion
	
	variable i
	for(i=0;i<num_cmp;i+=1)
		molarmass_ion[i]=12.000*ele_num[i][0]+1.008*ele_num[i][1]+15.995*ele_num[i][2]+14.003*ele_num[i][3]+126.904*ele_num[i][4]
	endfor
End

//Calculate the molar mass of each compound
Function GetCmpMass(Ele_num)
	wave Ele_num
	variable num_cmp=dimsize(ele_num,0)
	make/o/d/n=(num_cmp) MolarMass_cmp
	
	variable i
	for(i=0;i<num_cmp;i+=1)
		molarmass_cmp[i]=12.000*ele_num[i][0]+1.008*ele_num[i][1]+15.995*ele_num[i][2]+14.003*ele_num[i][3]
	endfor
End

//Calculate the molar mass of each compound
Function GetIonMass_s(Ele_num)
	wave Ele_num
	variable num_cmp=dimsize(ele_num,0)
	make/o/d/n=(num_cmp) MolarMass_ion
	
	variable i
	for(i=0;i<num_cmp;i+=1)
		molarmass_ion[i]=12.000*ele_num[i][0]+1.008*ele_num[i][1]+15.995*ele_num[i][2]+14.003*ele_num[i][3]+126.904*ele_num[i][4]+31.972*ele_num[i][5]
	endfor
End

//Calculate the molar mass of each compound
Function GetCmpMass_s(Ele_num)
	wave Ele_num
	variable num_cmp=dimsize(ele_num,0)
	make/o/d/n=(num_cmp) MolarMass_cmp
	
	variable i
	for(i=0;i<num_cmp;i+=1)
		molarmass_cmp[i]=12.000*ele_num[i][0]+1.008*ele_num[i][1]+15.995*ele_num[i][2]+14.003*ele_num[i][3]+31.972*ele_num[i][5]
	endfor
End


//After loading the text file, waves are named "waveXX"
//By default, wave0 is time wave, wave1 is temperature wave and the following waves are mass/signal of a compound
//This function renames wave0 and wave1, then creats a matrix for all the compounds and deletes all the wave "waveXX"
Function RenameWaves(num_cmp)
	variable num_cmp
	wave wave0,wave1
	
	duplicate/o wave0 time_s
	duplicate/o wave1 temp_C
	
	variable temp_num=numpnts(temp_c)
	make/o/d/n=(temp_num,num_cmp) cmp_matrix
	variable i
	string wave_str
	for(i=0;i<num_cmp;i+=1)
		wave_str="wave"+num2str(i+2)
		if(waveexists($wave_str))
			wave tempwave=$wave_str
			cmp_matrix[][i]=tempwave[p]
		else
			cmp_matrix[][i]=nan
		endif
	endfor
	
	for(i=0;i<num_cmp+2;i+=1)
		wave_str="wave"+num2str(i)
		killwaves $wave_str
	endfor
	
End

Function SeperateWavesNewFolder(seppnt,numwaves)
	variable seppnt
	variable numwaves
	wave wave0
	
	string currentfolder=Getdatafolder(1)
	string folderstr1=currentfolder+"Signal"
	string folderstr2=currentfolder+"Mass"
	Newdatafolder/o $folderstr1
	Newdatafolder/o $folderstr2
	
	variable temp_num=numpnts(wave0)
	variable deletepnts=temp_num-seppnt
	variable i
	string wave_str
	for(i=0;i<numwaves;i+=1)
		wave_str="wave"+num2str(i)
		if(waveexists($wave_str))
			duplicate/o $(currentfolder+wave_str) $(folderstr1+":"+wave_str)
			duplicate/o $(currentfolder+wave_str) $(folderstr2+":"+wave_str)			
		endif
	endfor
	
	for(i=0;i<numwaves;i+=1)
		wave_str="wave"+num2str(i)
		killwaves $wave_str
	endfor
	
	Setdatafolder $(folderstr1)
	for(i=0;i<numwaves;i+=1)
		wave_str="wave"+num2str(i)
		if(waveexists($wave_str))
			Deletepoints seppnt,deletepnts,$wave_str		
		endif
	endfor
	
	Setdatafolder $(folderstr2)
	for(i=0;i<numwaves;i+=1)
		wave_str="wave"+num2str(i)
		if(waveexists($wave_str))
			Deletepoints 0,seppnt,$wave_str		
		endif
	endfor
	
End


//This is different from  function RenameWaves to deal with data from Sigi
Function RenameWaves_s()
	wave wave0,wave1
	
	duplicate/o wave0 time_s
	duplicate/o wave1 temp_C
	wave CmpN_origin
	
	variable num_cmp=numpnts(cmpN_origin)
	variable temp_num=numpnts(temp_c)
	make/o/d/n=(temp_num,num_cmp) cmp_matrix
	variable i,cmpN
	string wave_str
	for(i=0;i<num_cmp;i+=1)
		cmpN=cmpN_origin[i]
		wave_str="wave"+num2str(cmpN-1)
		if(waveexists($wave_str))
			wave tempwave=$wave_str
			cmp_matrix[][i]=tempwave[p]
		else
			cmp_matrix[][i]=nan
		endif
	endfor
	
	for(i=0;i<2218;i+=1)
		wave_str="wave"+num2str(i)
		killwaves $wave_str
	endfor
	
End


//This function finds the max value of every column and normalizes every column to the max value
Function NormalizetoMax_col(matrix_origin)
	wave matrix_origin
	
	variable i
	string newwave_str
	newwave_str=nameofwave(matrix_origin)+"_norm_col"
	variable rownum=dimsize(matrix_origin,0)
	variable colnum=dimsize(matrix_origin,1)
	make/o/d/n=(rownum,colnum) $newwave_str
	make/o/d/n=(rownum)/Free tempwave
	wave newwave=$newwave_Str
	
	for(i=0;i<colnum;i+=1)
		tempwave=matrix_origin[p][i]
		wavestats/q tempwave
		tempwave/=V_max
		newwave[][i]=tempwave[p]
	endfor
	
	string notestr
	notestr="Normalized to the max values of every column of the matrix ' "+nameofwave(matrix_origin)+" ' "
	Note/k newwave notestr
End

//This function finds the max value of every row and normalizes every column to the max value
Function NormalizetoMax_row(matrix_origin)
	wave matrix_origin
	
	variable i
	string newwave_str
	newwave_str=nameofwave(matrix_origin)+"_norm_row"
	variable rownum=dimsize(matrix_origin,0)
	variable colnum=dimsize(matrix_origin,1)
	make/o/d/n=(rownum,colnum) $newwave_str
	make/o/d/n=(colnum)/Free tempwave
	wave newwave=$newwave_Str
	
	for(i=0;i<rownum;i+=1)
		tempwave=matrix_origin[i][p]
		wavestats/q tempwave
		tempwave/=V_max
		newwave[i][]=tempwave[q]
	endfor
	
	string notestr
	notestr="Normalized to the max values of every row of the matrix ' "+nameofwave(matrix_origin)+" ' "
	Note/k newwave notestr
End

//Plots every column
Function plotcol(ymatrix,xwave)
	wave ymatrix,xwave
	variable len=dimsize(ymatrix,1)
	variable i
	for(i=0;i<len;i+=1)
		if(i==0)
			display ymatrix[][i] vs xwave
		else
			appendtograph ymatrix[][i] vs xwave
		endif
	endfor
End

//Plots every row
Function plotrow(ymatrix,xwave)
	wave ymatrix,xwave
	variable len=dimsize(ymatrix,0)
	variable i
	for(i=0;i<len;i+=1)
		if(i==0)
			display ymatrix[i][] vs xwave
		else
			appendtograph ymatrix[i][] vs xwave
		endif
	endfor
End

//Plots every column
Function plotcol_filter(ymatrix,xwave,filterwave)
	wave ymatrix,xwave,filterwave
	variable len=dimsize(ymatrix,1)
	if(len!=numpnts(filterwave))
		print "Pls make sure that length of the filterwave is the same as the column number of ymatrix"
		abort
	endif
	variable i,count
	count=0
	for(i=0;i<len;i+=1)
		if(numtype(filterwave[i])==0)
			if(count==0)
				display ymatrix[][i] vs xwave
			else
				appendtograph ymatrix[][i] vs xwave
			endif
			count+=1
		endif
	endfor
End


//By default, this function integrates mass of each compound (represented by each column) over time
//ymatrix should be a 2D matrix containing mass/time of each compound at every time step
//xwave should be a 1D wave containing the time wave
//!!!!This function converts all the nan and negative values in the ymatrix to 0 to facilitate integration
Function IntegrateEachCol1(ymatrix,xwave)
	wave ymatrix, xwave
	
	variable rownum=dimsize(ymatrix,0)
	variable colnum=dimsize(ymatrix,1)
	variable wavelen=numpnts(xwave)
	if(rownum!=wavelen)
		print "make sure that the matrix has the row number same as num pnts of the x wave"
		abort
	endif
	
	string newxwave_str
	newxwave_str=nameofwave(xwave)+"_im"
	string newymatrix_str
	newymatrix_str=nameofwave(ymatrix)+"_INT"
	make/o/d/n=(wavelen+1) $newxwave_str
	wave newxwave=$newxwave_Str
	variable i,j
	//make the new xwave to facilitate integration
	if(xwave[0]==0)
		for(i=0;i<wavelen;i+=1)
			newxwave[i]=xwave[i]
		endfor
		newxwave[wavelen]=2*xwave[wavelen-1]-xwave[wavelen-2]
	else
		for(i=0;i<wavelen;i+=1)
			newxwave[i+1]=xwave[i]
		endfor	
		newxwave[0]=0
	endif
	
	make/o/d/n=(rownum)/Free tempmass,tempMFR
	make/o/d/n=(colnum) Cmp_mass
	make/o/d/n=(rownum,colnum) $newymatrix_str
	make/o/d/n=(rownum,colnum) Cmp_MFR
	wave newymatrix=$newymatrix_str
	variable lastvalue
	
	for(i=0;i<colnum;i+=1)
		tempmass=ymatrix[p][i]
		wavestats/q tempmass
		if(V_numNaNs!=0)
			for(j=0;j<rownum;j+=1)
				if(numtype(tempmass[j])==2)
					tempmass[j]=0
				elseif(tempmass[j]<0)
					tempmass[j]=0
				endif
			endfor
		else
			for(j=0;j<rownum;j+=1)
				if(tempmass[j]<0)
					tempmass[j]=0
				endif
			endfor
		endif
		Integrate tempmass/X=newxwave/D=tempmass_INT
		newymatrix[][i]=tempmass_int[p]
		lastvalue=tempmass_int[rownum-1]
		Cmp_mass[i]=lastvalue
		tempMFR=1-tempmass_int/lastvalue
		Cmp_MFR[][i]=tempMFR[p]
	endfor	
	
	string notestr
	notestr="This matrix stores integration of every column of ' "+nameofwave(ymatrix)+" ' over time wave ' "+nameofwave(newxwave)+" ' "
	Note/k newymatrix notestr
	notestr="This wave stores the integrated mass of every compound, the last row of matrix ' "+newymatrix_str+" ' "
	Note/k Cmp_mass notestr
	notestr="This matrix stores the MFR of every compound calculated from ' "+newymatrix_str+" ', using formula 1-Int/Int_max "
	Note/k Cmp_MFR notestr
	
End

//By default, this function adds up mass of all the compounds (represented by each column) at each time step
//ymatrix should be a 2D matrix containing mass of each compound at every time step
//!!!!This function converts all the nan and negative values in the ymatrix to 0 to facilitate integration
Function SumEachRow_noNan(ymatrix)
	wave ymatrix
	wave ele_num
	
	variable rownum=dimsize(ymatrix,0)
	variable colnum=dimsize(ymatrix,1)
	
	string newwave_str
	string parent_str
	parent_str=getdatafolder(0)
	if(stringmatch(parent_str,"Mass"))
		newwave_str="Totalmass_time"
	elseif(stringmatch(parent_Str,"Signal"))
		newwave_str="Totalsignal_time"
	else
		print "Something is wrong..."
		newwave_str="Total_time"
	endif

	make/o/d/n=(rownum) $newwave_str
	
	variable i,j
	
	make/o/d/n=(colnum)/Free tempmass,tempMFR
	make/o/d/n=(rownum) $newwave_str
	wave newwave=$newwave_str
	
	for(i=0;i<rownum;i+=1)
		tempmass=ymatrix[i][p]
		for(j=0;j<colnum;j+=1)
			if(numtype(tempmass[j])==2)
				tempmass[j]=0
			elseif(ele_num[j][0]==0 || ele_num[j][1]==0 || ele_num[j][2]==0)
				tempmass[j]=0
			//elseif(ele_num[j][5]!=0)
				//tempmass[j]=0
			//elseif(tempmass[j]<0)
				//tempmass[j]=0
			endif
		endfor
		wavestats/q tempmass
		newwave[i]=V_sum
	endfor	
	
	string notestr
	notestr="This matrix stores sum of every row of ' "+nameofwave(ymatrix)+" ' as a function of time "
	Note/k newwave notestr
	
End


//This function filters out ions that contribute less than xx% of the total mass
//Make sure that the originwave contains values in fraction, not percentage
Function FilterbyFraction(originwave,fraction,targetwave)
	wave originwave
	variable fraction  
	wave targetwave
	
	if(numpnts(originwave)!=numpnts(targetwave))
		print "Pls make sure that originwave and targetwave have the same length"
		abort
	endif
	
	string newwave_str,newtarget_str
	newwave_str=nameofwave(originwave)+"_filter"
	newtarget_str=nameofwave(targetwave)+"_filter"
	duplicate/o originwave $newwave_str
	duplicate/o targetwave $newtarget_str
	wave newwave=$newwave_str
	wave newtarget=$newtarget_str
	variable len=numpnts(originwave)
	variable i,count
	count=0
	
	for(i=0;i<len;i+=1)
		if(originwave[i]>=fraction)
			newwave[i]=originwave[i]
			newtarget[i]=targetwave[i]
			count+=1
		else
			newwave[i]=nan
			newtarget[i]=nan
		endif
	endfor
	
	string notestr
	notestr="This wave filters out the pnts in ' "+nameofwave(originwave)+" ' that have values below "+num2str(fraction)+"; Theere are "+num2str(count)+" pnts left after filtration"
	Note/k newwave notestr
	
End

Function CalculateHC()
	wave ele_num
	variable cmp_num=dimsize(ele_num,0)
	make/o/d/n=(cmp_num) HC_ratio,OC_ratio,NC_ratio,OSc
	
	variable i
	for(i=0;i<cmp_num;i+=1)
		if(ele_num[i][0]!=0)
			HC_ratio[i]=ele_num[i][1]/ele_num[i][0]
			OC_ratio[i]=ele_num[i][2]/ele_num[i][0]
			NC_ratio[i]=ele_num[i][3]/ele_num[i][0]
			OSc[i]=2*OC_ratio[i]-HC_ratio[i]-5*NC_ratio[i]
		else
			HC_ratio[i]=nan
			OC_ratio[i]=nan
			NC_ratio[i]=nan
			OSc[i]=nan
		endif
	endfor
	
End

