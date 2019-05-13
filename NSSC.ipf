#pragma rtGlobals=3		// Use modern global access method and strict wave access.

// Functions to cluster data using the noise-sorted scanning clustering algorithm, developed by Ziyue Li and Chris Cappa
// Initial code written by Ziyue Li (2018-19)
// Some functions are specific to FIGAERO-CIMS data. Future aim to generalize. 

// Version updates
// v1.0 - initiall functions written


//******************************************************************************
//!!! Some functions are paired with suffix "_np". In this case, functions w/o _np have pop-up menus
//	  while functions w/ _np don't have pop-up menus 


//******************************************************************************
//This is a suggested order of using the functions in the menu to process and cluster the data set
//The following steps serve as a guide of how to MANUALLY process and cluster the data set
// Suggested Procedures

	
	//1. Examine data quality
	FindRamping() //Determine the ramping period based on T vs. time
	//Find the weird ions from two aspects
	FindWeirdIons_subrange()  //Find the thermograms with an extremely negative valley 
	FindWeirdIons_Negsum()  //Find the thermograms that have a negative sum of signals
	CombineIndx()	//Combine the above two results
	GraphIndx_col() //Plot thermograms in one graph according to selected indexwave
	GraphIndx_massSpec() //This also helps find determine which ions to remove
	RemoveWeirdIons() //Remove weird ions from the matrix according the selected indexwave
		//Run the removeWeidIons function several times for different waves you would like to refine, such as cmp_matrix, cmp_mass and Names_origin

	//2. Data Pre-processing
	//2.1. Create waves of basic propertied of the ions
	Createwave_BasicP() //Create basic waves regarding the data set, including mass spectrum, bulk thermogram, m/z wave
	//2.2. Normalization and down-weighting
	MovingAverageMatrix_Col() //Smooth the thermograms for every ion
	MovingAverageWave() //Run twice to smooth both the desorption time and temperature,  use the same average points as the matrix
	Normalize2SmoothedMax() //Normalize thermograms to the smoothed maximum
	AverageMatrix_Subrange() //Average only the soaking period of the thermograms for down-weighting
	AverageWave_Subrange() //Run twice to average the soaking period of the desorption time and temperture wave for down-weighting
	//2.3. Find the noise and mass threshold
	MovingAverageMatrix_Col() //Smooth the normalized thermograms for every ion
	CalculateNoise_Subrange() //Calculate the average difference between unsmoothed and smoothed thermograms during ramping period as noise
		//Note that the subrange is determined from the ramping period in smoothed thermogram
	FindNoiseThreshold() //Determine the threshold of noise to filter low-noise ions
	FindMassFraThreshold() //Determine the threshold of "significant" ions based on high mass contribution
	//2.4. Filter out high-noise ions
	FilterMatrix_Col() //Filter out thermograms of high-noise ions from the matrix
	Filter1DWave() //Filter out high-noise ions from waves such as mass contribution, noise level and compound index..
	FilterTextWave() //Filter out high-noise ions in the text wave of compound names
	Createwave_BasicP_ext() //Calculates properties such as # of C, H, O, N, OSc, H/C, O/C for filtered ions
	
	//3. Clustering
	//3.1. Find the optimal eps
	FindOptEps_NSSC() //Scaning through several eps for clustering to help find the optimal eps
	GraphOptEpsParas() //Graph the parameters as a function of eps to help find the optimal eps
	GraphTG0_fromMat() //Graph the TGs of clusters using a selected eps to examine and compare
	//3.2. Use the user-defined optimal eps to conduct the final clustering
	NSSC_w1M() //Use the optimal eps to cluster again
	SortClusters_2c() //Sort clusters by T_50mass
	CountClstSize() //Count the size (number of ions in each cluster) of clusters, including the unclustered
	CalculateClstMass() //Calculate the mass fractional contribution of clusters, including the unclustered
	ClusteronClusterNum_unWt() //Calculate the unweighted avg TG of clusters
	ClusteronClusterNum_Wt() //Calculate the weighted avg TG of clusters
	//3.3. Graph clustering results
	GraphClstTG2_OneGraph() //Plot the TGs (both wt and unwt avg and individual members) of clusters in one graph
	
	//4. Properties of clusters
	ClusterSpectraMatrix() //Make a matrix containing the spectrum of clusters
	FilterMatrix_Col() //Filter out noisy ions from the original (no normalization or down-weighting) TG matrix with absolute signals/mass
	SumThermogramClst() //Use the above matrix to make a matrix containing sum TGs of clusters
	AverageClstProperties() //Average properties such as H/C, Osc for clusters
	

//*******************************************************************************
//The following functions are a partly automated procedure to conduct the data processing and clusteirng

//-------------------------------------------------------------------------------
//STEP 1, exams the quality of data and removes weird ions
//This step involves manual operation
//By default, the current folder must contain four waves: Names_origin(N), time_s(M),temp_c(M),cmp_matrix(M * N)
//By default, the matrix cmp_matrix has thermograms of ions in columns
Function StandardProcess_Step1() 
	
	wave/t Names_origin
	wave temp_c, time_s, cmp_matrix
	
	//Make sure that current folder has the essential waves with correct dimensions
	
	if(waveexists(names_origin)==0 || waveexists(temp_c) ==0 || waveexists(time_s) == 0|| waveexists(cmp_matrix) == 0)
		DoAlert 0, "Please make sure that waves Names_origin, time_s,temp_c,cmp_matrix are in the current folder"
		abort
	endif
	if(numpnts(time_s)!=numpnts(temp_c))
		DoAlert 0,"Please make sure that the time and temperature waves have the same number of points"
		abort
	endif
	if(numpnts(time_s)!=dimsize(cmp_matrix,0))
		DoAlert 0,  "Please make sure that the number of rows in cmp_matrix equals the number of points of the time wave"
		abort
	endif
	if(numpnts(names_origin)!=dimsize(cmp_matrix,1))
		DoAlert 0, "Please make sure that the number of columns in cmp_matrix equals the number of points of the names_origin text wave"
		abort
	endif
	
	variable Pnt_EndRamp = FindRamping_np(temp_c, time_s)
	variable SN_multiplier = 3
	variable Pnt_avg = 100
	
	DoAlert 1,"Would you like to manually input the point where ramping period ends?"
	if(V_flag==1)
		prompt Pnt_endramp,"Enter the point"
		Doprompt "Enter values", Pnt_endramp
		
		if(V_flag)
			abort
		endif
	endif
	
	//Find weird ions and store the index of those ions in a 1D wave
	FindWeirdIons_subrange_np(cmp_matrix,SN_multiplier,Pnt_avg,Pnt_avg/2,Pnt_EndRamp+150)
	wave TooNeg_indx
	FindWeirdIons_Negsum_np(cmp_matrix)
	wave cmp_mass
	wave cmp_mass_Neg_indx
	CombineIndx_np(TooNeg_indx,cmp_mass_Neg_indx,"Combined_indx")
	wave Combined_indx
	
	//Graph those weird and/or negative ions
	variable npanels_col 
	GraphIndx_massSpec_np(cmp_mass,Combined_indx)
	npanels_col = ceil((numpnts(Combined_indx))^0.5)
	GraphIndx_col_np(cmp_matrix,time_s,Combined_indx,npanels_col)
	
	make/O/I/U/N=0 Neg2Remove_indx
	Note/k Neg2Remove_indx "This wave stores the index of ions to be removed from the original data set"
	edit Neg2Remove_indx
	
	DoAlert 0,"Now please enter the index of weird ions based on the graphs into wave 'Neg2Remove_indx' "
	
End

//-------------------------------------------------------------------------------
//STEP 2, Preprocessing of the data
//This step is completely automatic
Function StandardProcess_Step2() 

	DoAlert 1,"Have you completed Step 1? "
	
	if(V_flag==2)
		abort
	endif
	
	DoAlert 1,"Have you set the current folder to the folder containing the original (unrefined) data set? "
	
	if(V_flag==2)
		abort
	endif
	
	string Datafldr_org = GetDataFolder(1)
	NewDatafolder/O RefinedData
	string Datafldr_rf = Datafldr_org+"RefinedData:"
	
	wave/t Names_origin
	wave time_s, temp_c
	wave cmp_matrix
	wave Neg2Remove_indx
	
	RemoveWeirdIons_np(cmp_matrix,Neg2Remove_indx)
	RemoveWeirdIons_np(Names_origin,Neg2Remove_indx)
	
	//Copy essential waves to the folder "RefinedData"
	duplicate/o $"time_s" $(Datafldr_rf+"time_s")
	duplicate/o $"temp_c" $(Datafldr_rf+"temp_c")
	duplicate/o $(Datafldr_org+"cmp_matrix_rf"), $(Datafldr_rf+"cmp_matrix_rf")
	duplicate/o $(Datafldr_org+"Names_origin_rf"), $(Datafldr_rf+"Names_origin_rf")
	
	SetDatafolder Datafldr_rf
	
	wave cmp_matrix_rf
	wave/t Names_origin_rf
	
	variable Pnt_EndRamp = FindRamping_np(temp_c, time_s)
	
	DoAlert 1,"Would you like to manually input the point where ramping period ends?"
	if(V_flag==1)
		prompt Pnt_endramp,"Enter the point"
		Doprompt "Enter values", Pnt_endramp
		
		if(V_flag)
			abort
		endif
	endif
	
	variable Pnt_EndSoak = numpnts(time_s)-1
	
	//Create some basic waves such as bulk TG, mass spectrum, m/z wave
	Createwave_BasicP_np(cmp_matrix_rf, Names_origin_rf,"_rf")
	wave Molarmass_cmp_rf
	wave cmp_mass_rf
	wave cmp_mass_norm_rf
	wave cmp_indx_rf
	
	variable Pnt_smtavg = 35
	variable Pnt_avg = 10
	variable Pnt_EndRamp_smt = Pnt_EndRamp-floor(Pnt_smtavg/2)
	
	//Smoothing
	MovingAverageMatrix_Col_np(cmp_matrix_rf,Pnt_smtavg)
	MovingAverageWave_np(time_s, Pnt_smtavg)
	MovingAverageWave_np(temp_c, Pnt_smtavg)
	wave cmp_matrix_rf_Smt35 = $("cmp_matrix_rf_Smt"+num2istr(Pnt_smtavg))
	wave time_s_smt35 = $("time_s_Smt"+num2istr(Pnt_smtavg))
	wave temp_c_smt35 = $("temp_c_Smt"+num2istr(Pnt_smtavg))
	
	//Normalization
	Normalize2SmoothedMax_np(cmp_matrix_rf_smt35,cmp_matrix_rf)
	wave cmp_matrix_rf_Smtnorm
	
	//Down-weighting
	AverageMatrix_Subrange_np(cmp_matrix_rf_smtnorm,Pnt_avg,Pnt_EndRamp+1,Pnt_EndSoak,"_10a")
	AverageWave_Subrange_np(time_s,Pnt_avg,Pnt_EndRamp+1,Pnt_EndSoak,"_10a")
	AverageWave_Subrange_np(temp_c,Pnt_avg,Pnt_EndRamp+1,Pnt_EndSoak,"_10a")
	wave cmp_matrix_rf_smtnorm_10a
	wave time_s_10a, temp_c_10a
	
	//Calculate the noise levels
	MovingAverageMatrix_Col_np(cmp_matrix_rf_smtnorm,Pnt_smtavg)
	wave cmp_matrix_rf_Smtnorm_Smt35 = $("cmp_matrix_rf_Smtnorm_Smt"+num2istr(Pnt_smtavg))
	CalculateNoise_Subrange_np(cmp_matrix_rf_smtnorm,cmp_matrix_rf_smtnorm_smt35,temp_c,temp_c_smt35,0,Pnt_EndRamp_smt,"_ar35")
	wave AvgNoise_ar35
	
	//Determine the noise threshold
	variable Ionpercent = 5
	variable multiplier = 3
	FindNoiseThreshold_np(AvgNoise_ar35,cmp_mass_norm_rf,Ionpercent,multiplier)
	wave Filterwave_noise	
	
	//Filtering out high-noise ions from the waves
	 FilterMatrix_Col_np(cmp_matrix_rf_smtnorm_10a,filterwave_noise,"cmp_matrix_sn_10a_fil")
	 Filter1DWave_np(cmp_mass_rf,filterwave_noise,"cmp_mass_fil")
	 Filter1DWave_np(cmp_mass_norm_rf,filterwave_noise,"cmp_mass_norm_fil")
	 Filter1DWave_np(cmp_indx_rf,filterwave_noise,"cmp_indx_fil")
	 Filter1DWave_np(AvgNoise_ar35,filterwave_noise,"AvgNoise_ar35_fil")
	 FilterTextWave_np(Names_origin_rf,filterwave_noise,"Names_fil")
	 wave/t Names_fil
	 
	 Createwave_BasicP_ext_np(Names_fil,"_fil")
	
End


//-------------------------------------------------------------------------------
//STEP 3.1, Clustering of the data -- Finding the optimal eps
//This step requires some manual operation, including examining the parameters and plot clustering results
Function StandardProcess_Step3_1() 

	DoAlert 1,"Have you completed Step 2? "
	
	if(V_flag==2)
		abort
	endif
	
//Determine the mass contribution threshold
	wave cmp_mass_norm_rf
	variable cumMassfra = 0.8
	variable HM_threshold = FindMassFraThreshold_np(cmp_mass_norm_rf,cumMassfra)
	
	wave cmp_matrix_sn_10a_fil
	wave time_s, temp_c
	wave time_s_10a, temp_c_10a
	wave cmp_mass_norm_fil
	wave Avgnoise_ar35_fil
	
	variable Pnt_EndRamp = FindRamping_np(temp_c, time_s)
	
	//Scan through a range of epsilon to get parameters such as Nc_all, Nc_no1M, fm_unclustered, R_avginterED
	FindOptEps_NSSC_np(cmp_matrix_sn_10a_fil,1.0,0.1,30,Pnt_EndRamp,HM_threshold,cmp_mass_norm_fil,temp_c_10a,AvgNoise_ar35_fil,"")
	wave epsilon
	wave ungroupedMass
	wave NumofClusters_all
	wave NumofClusters_no1M
	wave Ratio_avginterSD_eps
	
	GraphOptEpsParas_np(NumofClusters_all,Numofclusters_no1M,ungroupedMass,Ratio_avginterSD_eps,epsilon)
	
	DoAlert 0,"Now please examine the parameters to determine the optimal eps to use"
	
End

//-------------------------------------------------------------------------------
//STEP 3.2, Clustering of the data -- Final clustering and sorting
//This step is completely automated
Function StandardProcess_Step3_2() 
	
	DoAlert 1,"Have you completed Step 3.1? "
	
	if(V_flag==2)
		abort
	endif
	
	variable eps = 3.0
	
	prompt eps, "Enter the optimal eps to use"
	Doprompt "Enter values", eps
	
	wave cmp_matrix_sn_10a_fil
	
	CalculateED_np(cmp_matrix_sn_10a_fil,"ED_fil")
	wave ED_fil
	
	wave AvgNoise_ar35_fil
	wave cmp_mass_norm_fil
	wave cmp_mass_norm_rf
	wave time_s, temp_c
	wave time_s_10a, temp_c_10a
	
	variable Pnt_EndRamp = FindRamping_np(temp_c, time_s)
	variable cumMassfra = 0.8
	variable HM_threshold = FindMassFraThreshold_np(cmp_mass_norm_rf,cumMassfra)
	variable presetNc = ceil(dimsize(cmp_matrix_sn_10a_fil,1)/3)
	
	NSSC_w1M_np(cmp_matrix_sn_10a_fil,ED_fil,AvgNoise_ar35_fil,cmp_mass_norm_fil,presetNc,eps,HM_threshold)
	wave clusterave_wt
	wave clusternum
	wave clusterseedindx
	
	SortClusters_2c_np(Clusterave_wt,Clusternum,Clusterseedindx,temp_c_10a,"sort")
	wave Clusternumsort
	
	//Calculate unweighted and weighted average TG of sorted clusters
	ClusteronClusterNum_unWt_np(Clusternumsort,cmp_matrix_sn_10a_fil,"_unwtsort")
	ClusteronClusterNum_Wt_np(Clusternumsort,cmp_matrix_sn_10a_fil,cmp_mass_norm_fil,"_wtsort")
	wave Clusterave_unwtsort, Clusterave_wtsort
	
	//Calculate the size and mass fractional contributions of sorted clusters
	CountClstSize_np(Clusternumsort,"Clustersizesort")
	CalculateClstMass_np(Clusternumsort,cmp_mass_norm_fil,"Clustermasssort")
	
	variable npanels_col=ceil((dimsize(clusterave_wtsort,1))^0.5)
	
	//Plots the unweighted and weighted avg TG of clusters, along with the TGs of individual members
	GraphClstTG2_OneGraph_np(clusternumsort, cmp_matrix_sn_10a_fil,clusterave_unwtsort,clusterave_wtsort, temp_c_10a, npanels_col)

End

//-------------------------------------------------------------------------------
//STEP 4, Average properties of each cluster
//This step is completely automated
//This step creates many graphs
Function StandardProcess_Step4()
	
	DoAlert 1,"Have you completed Step 3.2? "
	
	if(V_flag==2)
		abort
	endif
	
	wave cmp_mass_norm_fil
	wave cmp_mass_fil
	wave Molarmass_cmp_fil
	wave clusternumsort, clustersizesort, clustermasssort
	wave filterwave_noise
	wave cmp_matrix_rf
	wave time_s, temp_c
	wave time_s_10a, temp_c_10a
	wave Mass50Locsort
	wave clusterave_wtsort
	
	ClusterSpectraMatrix_np(cmp_mass_fil,clusternumsort,"_abssort")
	ClusterSpectraMatrix_np(cmp_mass_norm_fil,clusternumsort,"_normsort")
	wave ClusterSpectra_sort
	wave ClusterSpectra_normsort
	
	FilterMatrix_Col_np(cmp_matrix_rf,filterwave_noise,"cmp_matrix_fil")
	wave cmp_matrix_fil
	
	SumThermogramClst_np(cmp_matrix_fil,clusternumsort,"sort")
	wave SumTG_clstsort
	wave Bulk_TG_rf
	
	AverageClstProperties_np(Clusternumsort, Clustersizesort)
	
	Loc2T_np(Mass50Locsort,temp_c_10a,"Mass50Tsort")
	
	GraphMatrix_Col_np(clusterave_wtsort, temp_c_10a, "ClusterAve_wt")
	ModifyGraph width=250, height=150 
	ModifyGraph mirror=1,standoff=0
	Label left "Normalized mass"
	Label bottom "Desorption temperature, ¡C"
	SetAxis left -0.1,1.3
	ModifyGraph lsize=1.5
	GraphClstSpectrum_stacked_np(ClusterSpectra_normsort,molarmass_cmp_fil)
	GraphClstTG_stacked_np(SumTG_clstsort,Bulk_TG_rf,time_s)
	GraphClstSize_np(Clustersizesort)
	GraphClstMass_np(Clustermasssort)
	
End 

//-------------------------------------------------------------------------------
//-------------------------------------------------------------------------------
//STEP 1, Refine and filter waves based on existing indexwave and filterwave
//This step needs input of certain waves
Function StandardProcess1_Step1() 

	wave/t Names_origin
	wave temp_c, time_s, cmp_matrix
	
	//Make sure that current folder has the essential waves with correct dimensions
	
	if(waveexists(names_origin)==0 || waveexists(temp_c) ==0 || waveexists(time_s) == 0|| waveexists(cmp_matrix) == 0)
		DoAlert 0, "Please make sure that waves Names_origin, time_s,temp_c,cmp_matrix are in the current folder"
		abort
	endif
	if(numpnts(time_s)!=numpnts(temp_c))
		DoAlert 0,"Please make sure that the time and temperature waves have the same number of points"
		abort
	endif
	if(numpnts(time_s)!=dimsize(cmp_matrix,0))
		DoAlert 0,  "Please make sure that the number of rows in cmp_matrix equals the number of points of the time wave"
		abort
	endif
	if(numpnts(names_origin)!=dimsize(cmp_matrix,1))
		DoAlert 0, "Please make sure that the number of columns in cmp_matrix equals the number of points of the names_origin text wave"
		abort
	endif
	
	variable Pnt_EndRamp = FindRamping_np(temp_c, time_s)
	variable SN_multiplier = 3
	variable Pnt_avg = 100
	
	DoAlert 1,"Would you like to manually input the point where ramping period ends?"
	if(V_flag==1)
		prompt Pnt_endramp,"Enter the point"
		Doprompt "Enter values", Pnt_endramp
		
		if(V_flag)
			abort
		endif
	endif
	
	string Negindexstr="Neg2remove_indx"
	String filterwavestr="filterwave_noise"
	String Clusternumstr="Clusternumsort"
	
	prompt Negindexstr, "Select index wave of weird ions", popup WaveList("*",";","")
	prompt filterwavestr, "Select filter wave of noise", popup WaveList("*",";","")
	prompt Clusternumstr, "Select wave of clst# for ions", popup WaveList("*",";","")
	Doprompt "Select waves",Negindexstr, filterwavestr, clusternumstr
	
	if(V_flag)
		abort
	endif
	
	string Datafldr_org = GetDataFolder(1)
	NewDatafolder/O RefinedData
	string Datafldr_rf = Datafldr_org+"RefinedData:"
	
	wave Neg2Remove_indx=$Negindexstr
	wave filterwave_noise=$filterwavestr
	wave clusternum=$clusternumstr
	
	RemoveWeirdIons_np(cmp_matrix,Neg2Remove_indx)
	RemoveWeirdIons_np(Names_origin,Neg2Remove_indx)
	
	//Copy essential waves to the folder "RefinedData"
	duplicate/o $"time_s" $(Datafldr_rf+"time_s")
	duplicate/o $"temp_c" $(Datafldr_rf+"temp_c")
	duplicate/o $(Datafldr_org+"cmp_matrix_rf"), $(Datafldr_rf+"cmp_matrix_rf")
	duplicate/o $(Datafldr_org+"Names_origin_rf"), $(Datafldr_rf+"Names_origin_rf")
	duplicate/o $(Datafldr_org+filterwavestr), $(Datafldr_rf+filterwavestr)
	duplicate/o $(Datafldr_org+clusternumstr), $(Datafldr_rf+clusternumstr)
	
	SetDatafolder Datafldr_rf
	
	wave cmp_matrix_rf
	wave/t Names_origin_rf
	
	variable Pnt_EndSoak = numpnts(time_s)-1
	
	//Create some basic waves such as bulk TG, mass spectrum, m/z wave
	Createwave_BasicP_np(cmp_matrix_rf, Names_origin_rf,"_rf")
	wave Molarmass_cmp_rf
	wave cmp_mass_rf
	wave cmp_mass_norm_rf
	wave cmp_indx_rf
	
	variable Pnt_smtavg = 35
	Pnt_avg = 10
	variable Pnt_EndRamp_smt = Pnt_EndRamp-floor(Pnt_smtavg/2)
	
	//Smoothing
	MovingAverageMatrix_Col_np(cmp_matrix_rf,Pnt_smtavg)
	MovingAverageWave_np(time_s, Pnt_smtavg)
	MovingAverageWave_np(temp_c, Pnt_smtavg)
	wave cmp_matrix_rf_Smt35 = $("cmp_matrix_rf_Smt"+num2istr(Pnt_smtavg))
	wave time_s_smt35 = $("time_s_Smt"+num2istr(Pnt_smtavg))
	wave temp_c_smt35 = $("temp_c_Smt"+num2istr(Pnt_smtavg))
	
	//Normalization
	Normalize2SmoothedMax_np(cmp_matrix_rf_smt35,cmp_matrix_rf)
	wave cmp_matrix_rf_Smtnorm
	
	//Down-weighting
	AverageMatrix_Subrange_np(cmp_matrix_rf_smtnorm,Pnt_avg,Pnt_EndRamp+1,Pnt_EndSoak,"_10a")
	AverageWave_Subrange_np(time_s,Pnt_avg,Pnt_EndRamp+1,Pnt_EndSoak,"_10a")
	AverageWave_Subrange_np(temp_c,Pnt_avg,Pnt_EndRamp+1,Pnt_EndSoak,"_10a")
	wave cmp_matrix_rf_smtnorm_10a
	wave time_s_10a, temp_c_10a
	
	//Calculate the noise levels
	MovingAverageMatrix_Col_np(cmp_matrix_rf_smtnorm,Pnt_smtavg)
	wave cmp_matrix_rf_Smtnorm_Smt35 = $("cmp_matrix_rf_Smtnorm_Smt"+num2istr(Pnt_smtavg))
	CalculateNoise_Subrange_np(cmp_matrix_rf_smtnorm,cmp_matrix_rf_smtnorm_smt35,temp_c,temp_c_smt35,0,Pnt_EndRamp_smt,"_ar35")
	wave AvgNoise_ar35	
	
	//Filtering out high-noise ions from the waves
	 FilterMatrix_Col_np(cmp_matrix_rf_smtnorm_10a,filterwave_noise,"cmp_matrix_sn_10a_fil")
	 FilterMatrix_Col_np(cmp_matrix_rf,filterwave_noise,"cmp_matrix_fil")
	 Filter1DWave_np(cmp_mass_rf,filterwave_noise,"cmp_mass_fil")
	 Filter1DWave_np(cmp_mass_norm_rf,filterwave_noise,"cmp_mass_norm_fil")
	 Filter1DWave_np(cmp_indx_rf,filterwave_noise,"cmp_indx_fil")
	 Filter1DWave_np(AvgNoise_ar35,filterwave_noise,"AvgNoise_ar35_fil")
	 FilterTextWave_np(Names_origin_rf,filterwave_noise,"Names_fil")
	 wave/t Names_fil
	 
	 Createwave_BasicP_ext_np(Names_fil,"_fil")
	 
	 duplicate/o clusternum Clusternumsort_org
	
End

//-------------------------------------------------------------------------------
//STEP 2, Clustering based on existing clusternum wave
//This step is completely automated
Function StandardProcess1_Step2() 
	
	DoAlert 1,"Have you completed Step 1? "
	
	if(V_flag==2)
		abort
	endif
	
	wave Clusternumsort_org
	
	wave cmp_matrix_sn_10a_fil
	wave temp_c_10a, time_s_10a
	wave cmp_mass_norm_fil
	
	//Calculate unweighted and weighted average TG of sorted clusters
	ClusteronClusterNum_unWt_np(Clusternumsort_org,cmp_matrix_sn_10a_fil,"_unwtsort")
	ClusteronClusterNum_Wt_np(Clusternumsort_org,cmp_matrix_sn_10a_fil,cmp_mass_norm_fil,"_wtsort")
	wave Clusterave_unwtsort, Clusterave_wtsort
	
	//Calculate the size and mass fractional contributions of sorted clusters
	CountClstSize_np(Clusternumsort_org,"Clustersizesort")
	CalculateClstMass_np(Clusternumsort_org,cmp_mass_norm_fil,"Clustermasssort")
	
	variable npanels_col=ceil((dimsize(clusterave_wtsort,1))^0.5)
	
	//Plots the unweighted and weighted avg TG of clusters, along with the TGs of individual members
	GraphClstTG2_OneGraph_np(clusternumsort_org, cmp_matrix_sn_10a_fil,clusterave_unwtsort,clusterave_wtsort, temp_c_10a, npanels_col)

	wave cmp_matrix_fil
	
	SumThermogramClst_np(cmp_matrix_fil,clusternumsort_org,"sort")

End


//*******************************************************************************
menu "Clustering"
	Submenu "Individual Clustering"
		"Step 1: Data quality",StandardProcess_Step1()
		"Step 2: Pre-processing",StandardProcess_Step2()
		"Step 3.1: Optimal eps",StandardProcess_Step3_1()
		"Step 3.2: Clustering",StandardProcess_Step3_2()
		"Step 4: Cluster properties",StandardProcess_Step4()
	end
	Submenu "Sequent Clustering"
		"Step 1: Pre-processing",StandardProcess1_Step1()
		"Step 2: Clustering", StandardProcess1_Step2()
	end
	Submenu "Desorption Period"
		"Determine ramping period", FindRamping()
	end	
	Submenu "Weird Ions"
		"Identify weird ions", FindWeirdIons_subrange()
		"Identify Negsum ions", FindWeirdIons_Negsum()
		"Removal", RemoveWeirdIons()
		"Graph TG", GraphIndx_col()
		"Graph MassSpectrum", GraphIndx_massSpec()
	end
	Submenu "Pre-Processing"
		"BasicProperties", Createwave_BasicP()
		"Smoothing matrix", MovingAverageMatrix_Col()
		"Smoothing wave", MovingAverageWave()
		"Average part of the matrix",  AverageMatrix_Subrange()
		"Average part of the wave",  AverageWave_Subrange()
		"Normalization", Normalize2SmoothedMax()
		"Calculate noise levels", CalculateNoise_Subrange()
		"Find noise threshold", FindNoiseThreshold()
		"Find massfraction threshold", FindMassFraThreshold()
		"Make filtered matrix", FilterMatrix_Col()
		"Make filtered wave", Filter1DWave()
	end
	Submenu "Clustering"
		"Calculate Euclidean distance", CalculateED()
		"NSSC without one-member clusters", NSSC()
		"NSSC with one-member clusters", NSSC_w1M()
		"Find the optimal eps", FindOptEps_NSSC()
		"Sorting of clusters", SortClusters_2c()
		"Count the size of clusters", CountClstsize()
		"Calculate the mass fraction of clusters", CalculateClstMass()
	end
	Submenu "Clst properties"
		"Separate mass spectrum of clusters", ClusterSpectraMatrix()
		"Sum Thermograms of clusters", SumThermogramClst()
		"Avg thermograms of clusters", ClusteronClusterNum_unWt()
		"Weighted avg thermograms of clusters", ClusteronClusterNum_Wt()
		"Calculate the avg properties of clusters", AverageClstProperties()
		"Calculate the weighted avg of selected property of clusters",AverageClst_specific_wt()
	end
	Submenu "Graphing"
		"Plot parameters vs eps, post FindOptEps", GraphOptEpsParas()
		"Plot TGs of clusters, post FindOptEps", GraphTG0_fromMat()
		"Plot TGs of clusters, avg and individual", GraphClstTG2_oneGraph()
		"Plot stacked TG of clusters and bulk TG", GraphClstTG_stacked()
		"Plot stacked spectra of clusters", GraphClstSpectrum_stacked()
	end
end

//*******************************************************************************


//This function aims to define the ramping period from soaking period based on the change of temperature with time
Function FindRamping()
	String tempwavestr="temp_c"
	String timewavestr="time_s"
	
	prompt tempwavestr, "Select temperature wave", popup WaveList("*",";","")
	prompt timewavestr, "Select time wave", popup WaveList("*",";","")
	Doprompt "Select waves",tempwavestr,timewavestr
	
	if(V_flag)
		abort
	endif
	
	wave tempwave=$tempwavestr
	wave timewave=$timewavestr
	
	if(numpnts(tempwave)!=numpnts(timewave))
		print "Please make sure the temperature and time wave have the same number of points"
		abort
	endif
	
	variable npnts=numpnts(timewave)
	
	variable i,j
	
	make/o/d/n=(npnts-1)/Free Slope2
	Slope2=0
	
	//Calculate the slopes between every two points
	for(i=0;i<npnts-1;i+=1)
		Slope2[i]=(tempwave[i+1]-tempwave[i])/(timewave[i+1]-timewave[i])
	endfor
	
	variable Avgnpnt=9  //the number of points of slopes to be averaged
	variable npnts2=npnts-Avgnpnt+1 
	
	make/o/d/n=(npnts2)/Free AvgSlope
	AvgSlope=0
	
	for(i=0;i<npnts2-1;i+=1)
		for(j=0;j<avgnpnt;j+=1)
			AvgSlope[i]+=Slope2[i+j]
		endfor
		AvgSlope[i]/=Avgnpnt
	endfor
	
	make/o/d/n=(100)/Free wave2avg_soak
	wave2avg_soak=0
	
	for(i=0;i<100;i+=1)
		wave2avg_soak[i]=Avgslope[npnts2-1-i]
	endfor
	wavestats/q wave2avg_soak
	
	variable AvgT_soak=V_avg
	variable SdvT_soak=V_sdev
	
	variable Threshold=AvgT_soak+5*SdvT_soak
	
	//Find backwards the point where the average slope suddenly increases
	//Because for both ramping and soaking periods, the slopes are rather flat
	//However the slope of ramping period is significantly higher than soaking period
	do
	
		Findlevel/EDGE=2/Q/P/R=[npnts2-1,0] AvgSlope,threshold
		variable EndofRamp
		EndofRamp=round(V_LevelX)+(Avgnpnt-1)/2
	
		make/o/d/n=(50)/Free wave2avg_soak1
		wave2avg_soak1=0
	
		for(i=0;i<50;i+=1)
			wave2avg_soak1[i]=Avgslope[EndofRamp-i]
		endfor
		wavestats/q wave2avg_soak1
	
		variable AvgT_temp = V_avg
		
		DeletePoints (EndofRamp-Avgnpnt),npnts, AvgSlope
	
	while(AvgT_temp < Threshold)
	
	print "The point where ramping period ends is:"
	print EndofRamp
	
End

Function FindRamping_np(tempwave, timewave)
	
	wave tempwave
	wave timewave
	
	if(numpnts(tempwave)!=numpnts(timewave))
		print "Please make sure the temperature and time wave have the same number of points"
		abort
	endif
	
	variable npnts=numpnts(timewave)
	
	variable i,j
	
	make/o/d/n=(npnts-1)/Free Slope2
	Slope2=0
	
	//Calculate the slopes between every two points
	for(i=0;i<npnts-1;i+=1)
		Slope2[i]=(tempwave[i+1]-tempwave[i])/(timewave[i+1]-timewave[i])
	endfor
	
	variable Avgnpnt=9  //the number of points of slopes to be averaged
	variable npnts2=npnts-Avgnpnt+1 
	
	make/o/d/n=(npnts2)/Free AvgSlope
	AvgSlope=0
	
	for(i=0;i<npnts2-1;i+=1)
		for(j=0;j<avgnpnt;j+=1)
			AvgSlope[i]+=Slope2[i+j]
		endfor
		AvgSlope[i]/=Avgnpnt
	endfor
	
	make/o/d/n=(100)/Free wave2avg_soak
	wave2avg_soak=0
	
	for(i=0;i<100;i+=1)
		wave2avg_soak[i]=Avgslope[npnts2-1-i]
	endfor
	wavestats/q wave2avg_soak
	
	variable AvgT_soak=V_avg
	variable SdvT_soak=V_sdev
	
	variable Threshold=AvgT_soak+5*SdvT_soak
	
	//print AvgT_soak
	//print SdvT_soak
	//print Threshold
	//abort
	
	//Find backwards the point where the average slope suddenly increases
	//Because for both ramping and soaking periods, the slopes are rather flat
	//However the slope of ramping period is significantly higher than soaking period
	do
	
		Findlevel/EDGE=2/Q/P/R=[npnts2-1,0] AvgSlope,threshold
		variable EndofRamp
		EndofRamp=round(V_LevelX)+(Avgnpnt-1)/2
	
		make/o/d/n=(50)/Free wave2avg_soak1
		wave2avg_soak1=0
	
		for(i=0;i<50;i+=1)
			wave2avg_soak1[i]=Avgslope[EndofRamp-i]
		endfor
		wavestats/q wave2avg_soak1
	
		variable AvgT_temp = V_avg
		
		DeletePoints (EndofRamp-Avgnpnt),npnts, AvgSlope
	
	while(AvgT_temp < Threshold)
	
	//print "The point where ramping period ends is:"
	return EndofRamp
		
	
End


//Weird ions are defined as ions with thermograms that have a valley with negative values beyond x times the normal standard deviation
Function FindWeirdIons()
	
	String ymatrixstr="cmp_matrix"
	variable multiplier=3
	variable npntstd = 120
	
	prompt ymatrixstr, "Select signal matrix", popup WaveList("*",";","")
	prompt multiplier, "Multiplier of sdev"
	prompt npntstd, "num of pnts to be averaged"
	Doprompt "Select waves",ymatrixstr,multiplier, npntstd
	
	if(V_flag)
		abort
	endif
	
	wave ymatrix=$ymatrixstr //each column of the matrix is the thermogram of a different ion, must use absolute signal/mass instead of normalized
	variable npntabs=npntstd //numpber of points used to calculate the average of the smallest absolute signals
	
	variable nrows=dimsize(ymatrix,0)
	variable ncols=dimsize(ymatrix,1)
	
	variable i,j
	variable minloc
	
	make/o/d/n=(ncols) NormalStd //Level of uncertainty/noise
	make/o/d/n=(npntstd)/Free tempwave
	make/o/d/n=(ncols) ValleyAvg //Level of the negative values of the valley
	make/o/d/n=(ncols) ValleyLoc //Location of the minimum value
	make/o/d/n=(npntabs)/Free tempwave1
	make/o/d/n=(nrows-npntabs)/Free tempwave2
	
	
	//Calculate the standard deviation
	//We use the points towards the end of soaking period because it is usually flat there
	for(i=0;i<ncols;i+=1)
		for(j=0;j<npntstd;j+=1)
			tempwave[j]=ymatrix[nrows-10-1-j][i]
		endfor
		wavestats/q tempwave
		Normalstd[i]=V_sdev
	endfor
	
	//Calculate the average of absolute signalsof the valley
	for(i=0;i<ncols;i+=1)
		for(j=0;j<nrows-npntabs;j+=1)
			tempwave2[j]=ymatrix[npntabs*0.5+j][i]
		endfor
		wavestats/q tempwave2
		minloc=V_minloc
		Valleyloc[i]=minloc+npntabs*0.5
		for(j=0;j<npntabs;j+=1)
			tempwave1[j]=ymatrix[minloc+j][i]
		endfor
		wavestats/q tempwave1
		ValleyAvg[i]=V_avg
	endfor
	
	//variable threshold=multiplier*
	
	extract/o/indx normalstd,TooNeg_indx,Valleyavg<multiplier*Normalstd*(-1)
	
	Note/k Normalstd "This wave stores the sdev of last 100 points for every column in matrix "+ymatrixstr	
	Note/k ValleyLoc "This wave stores the locations (point) where the minimum values of each column of matrix "+ymatrixstr+" are"
	Note/k ValleyAvg "This wave stores the average value of the valley (+/- 50 points of the minimum) of each column of matrix "+ymatrixstr	
	Note/k TooNeg_indx "This wave stores the index of the weird ions"
	
End

Function FindWeirdIons_np(ymatrix, multiplier, npntstd)
	
	Wave ymatrix
	variable multiplier
	variable npntstd 
	
	string ymatrixstr = nameofwave(ymatrix)

	variable npntabs=npntstd //numpber of points used to calculate the average of the smallest absolute signals
	
	variable nrows=dimsize(ymatrix,0)
	variable ncols=dimsize(ymatrix,1)
	
	variable i,j
	variable minloc
	
	make/o/d/n=(ncols) NormalStd //Level of uncertainty/noise
	make/o/d/n=(npntstd)/Free tempwave
	make/o/d/n=(ncols) ValleyAvg //Level of the negative values of the valley
	make/o/d/n=(ncols) ValleyLoc //Location of the minimum value
	make/o/d/n=(npntabs)/Free tempwave1
	make/o/d/n=(nrows-npntabs)/Free tempwave2
	
	
	//Calculate the standard deviation
	//We use the points towards the end of soaking period because it is usually flat there
	for(i=0;i<ncols;i+=1)
		for(j=0;j<npntstd;j+=1)
			tempwave[j]=ymatrix[nrows-10-1-j][i]
		endfor
		wavestats/q tempwave
		Normalstd[i]=V_sdev
	endfor
	
	//Calculate the average of absolute signalsof the valley
	for(i=0;i<ncols;i+=1)
		for(j=0;j<nrows-npntabs;j+=1)
			tempwave2[j]=ymatrix[npntabs*0.5+j][i]
		endfor
		wavestats/q tempwave2
		minloc=V_minloc
		Valleyloc[i]=minloc+npntabs*0.5
		for(j=0;j<npntabs;j+=1)
			tempwave1[j]=ymatrix[minloc+j][i]
		endfor
		wavestats/q tempwave1
		ValleyAvg[i]=V_avg
	endfor
	
	//variable threshold=multiplier*
	
	extract/o/indx normalstd,TooNeg_indx,Valleyavg<multiplier*Normalstd*(-1)
	
	Note/k Normalstd "This wave stores the sdev of last 100 points for every column in matrix "+ymatrixstr	
	Note/k ValleyLoc "This wave stores the locations (point) where the minimum values of each column of matrix "+ymatrixstr+" are"
	Note/k ValleyAvg "This wave stores the average value of the valley (+/- 50 points of the minimum) of each column of matrix "+ymatrixstr	
	Note/k TooNeg_indx "This wave stores the index of the weird ions"
	
End


//Weird ions are defined as ions with thermograms that have a valley with negative values beyond x times the normal standard deviation
//Please see function FindWeirdIons if you want to use all the points for Min search
Function FindWeirdIons_subrange()
	
	String ymatrixstr="cmp_matrix"
	variable multiplier = 3
	variable npntstd = 120
	variable startpnt = npntstd/2
	variable endpnt = 400
	
	prompt ymatrixstr, "Select signal matrix", popup WaveList("*",";","")
	prompt multiplier, "Multiplier of sdev"
	prompt npntstd, "Num of pnts to be averaged"
	prompt startpnt, "Start pnt for Min search"
	prompt endpnt, "End pnt for Min search"
	Doprompt "Select waves",ymatrixstr,multiplier, npntstd, startpnt, endpnt
	
	if(V_flag)
		abort
	endif
	
	wave ymatrix=$ymatrixstr //each column of the matrix is the thermogram of a different ion, must use absolute signal/mass instead of normalized
	variable npntabs=npntstd-40 //numpber of points used to calculate the average of the smallest absolute signals
	
	variable nrows=dimsize(ymatrix,0)
	variable ncols=dimsize(ymatrix,1)
	
	variable i,j
	variable minloc
	
	make/o/d/n=(ncols) NormalStd //Level of uncertainty/noise
	make/o/d/n=(npntstd)/Free tempwave
	make/o/d/n=(ncols) ValleyAvg //Level of the negative values of the valley
	make/o/d/n=(ncols) ValleyLoc //Location of the minimum value
	make/o/d/n=(npntabs)/Free tempwave1
	make/o/d/n=(endpnt-startpnt+1)/Free tempwave2
	
	
	//Calculate the standard deviation
	//We use the points towards the end of soaking period because it is usually flat there
	for(i=0;i<ncols;i+=1)
		for(j=0;j<npntstd;j+=1)
			tempwave[j]=ymatrix[nrows-10-1-j][i]
		endfor
		wavestats/q tempwave
		Normalstd[i]=V_sdev
	endfor
	
	//Calculate the average of absolute signalsof the valley
	for(i=0;i<ncols;i+=1)
		for(j=0;j<endpnt-startpnt+1;j+=1)
			tempwave2[j]=ymatrix[startpnt+j][i]
		endfor
		wavestats/q tempwave2
		minloc=V_minloc+startpnt
		Valleyloc[i]=minloc
		for(j=0;j<npntabs;j+=1)
			tempwave1[j]=ymatrix[minloc-npntabs/2+j][i]
		endfor
		wavestats/q tempwave1
		ValleyAvg[i]=V_avg
	endfor
	
	//variable threshold=multiplier*
	
	extract/o/indx normalstd,TooNeg_indx,Valleyavg<multiplier*Normalstd*(-1)
	
	Note/k Normalstd "This wave stores the sdev of last 100 points for every column in matrix "+ymatrixstr	
	Note/k ValleyLoc "This wave stores the locations (point) where the minimum values of each column of matrix "+ymatrixstr+" are between row "+num2str(startpnt)+" and row "+num2str(endpnt) 
	Note/k ValleyAvg "This wave stores the average value of the valley (+/- 50 points of the minimum) of each column of matrix "+ymatrixstr	
	Note/k TooNeg_indx "This wave stores the index of the weird ions"
	
End

Function FindWeirdIons_subrange_np(ymatrix,multiplier,npntstd,startpnt,endpnt)
	
	wave ymatrix
	variable multiplier
	variable npntstd
	variable startpnt
	variable endpnt
	
	string ymatrixstr=nameofwave(ymatrix) //each column of the matrix is the thermogram of a different ion, must use absolute signal/mass instead of normalized
	variable npntabs=npntstd-40 //numpber of points used to calculate the average of the smallest absolute signals
	
	variable nrows=dimsize(ymatrix,0)
	variable ncols=dimsize(ymatrix,1)
	
	variable i,j
	variable minloc
	
	make/o/d/n=(ncols) NormalStd //Level of uncertainty/noise
	make/o/d/n=(npntstd)/Free tempwave
	make/o/d/n=(ncols) ValleyAvg //Level of the negative values of the valley
	make/o/d/n=(ncols) ValleyLoc //Location of the minimum value
	make/o/d/n=(npntabs)/Free tempwave1
	make/o/d/n=(endpnt-startpnt+1)/Free tempwave2
	
	
	//Calculate the standard deviation
	//We use the points towards the end of soaking period because it is usually flat there
	for(i=0;i<ncols;i+=1)
		for(j=0;j<npntstd;j+=1)
			tempwave[j]=ymatrix[nrows-10-1-j][i]
		endfor
		wavestats/q tempwave
		Normalstd[i]=V_sdev
	endfor
	
	//Calculate the average of absolute signalsof the valley
	for(i=0;i<ncols;i+=1)
		for(j=0;j<endpnt-startpnt+1;j+=1)
			tempwave2[j]=ymatrix[startpnt+j][i]
		endfor
		wavestats/q tempwave2
		minloc=V_minloc+startpnt
		Valleyloc[i]=minloc
		for(j=0;j<npntabs;j+=1)
			tempwave1[j]=ymatrix[minloc-npntabs/2+j][i]
		endfor
		wavestats/q tempwave1
		ValleyAvg[i]=V_avg
	endfor
	
	//variable threshold=multiplier*
	
	extract/o/indx normalstd,TooNeg_indx,Valleyavg<multiplier*Normalstd*(-1)
	
	Note/k Normalstd "This wave stores the sdev of last 100 points for every column in matrix "+ymatrixstr	
	Note/k ValleyLoc "This wave stores the locations (point) where the minimum values of each column of matrix "+ymatrixstr+" are between row "+num2str(startpnt)+" and row "+num2str(endpnt) 
	Note/k ValleyAvg "This wave stores the average value of the valley (+/- 50 points of the minimum) of each column of matrix "+ymatrixstr	
	Note/k TooNeg_indx "This wave stores the index of the weird ions"
	
End


//This function finds the ions that have a negative summed-up signals/mass
Function FindWeirdIons_Negsum()

	String ymatrixstr="cmp_matrix"
	
	prompt ymatrixstr, "Select signal matrix", popup WaveList("*",";","")
	Doprompt "Select waves",ymatrixstr
	
	if(V_flag)
		abort
	endif
	
	wave ymatrix=$ymatrixstr 
	
	matrixop/o cmp_mass = sumcols(ymatrix)^t
	
	extract/o cmp_mass, cmp_mass_Neg, cmp_mass<0
	extract/o/indx cmp_mass, cmp_mass_Neg_indx, cmp_mass<0
	
	Note/k cmp_mass "This wave stores the sum of each column of matrix "+ymatrixstr
	Note/k cmp_mass_Neg "This wave contains the negative values in wave cmp_mass"
	Note/k cmp_mass_Neg_indx "This wave contains the index number of negative values in wave cmp_mass"
End

Function FindWeirdIons_Negsum_np(ymatrix)
	wave ymatrix

	String ymatrixstr=nameofwave(ymatrix)
	
	matrixop/o cmp_mass = sumcols(ymatrix)^t
	
	extract/o cmp_mass, cmp_mass_Neg, cmp_mass<0
	extract/o/indx cmp_mass, cmp_mass_Neg_indx, cmp_mass<0
	
	Note/k cmp_mass "This wave stores the sum of each column of matrix "+ymatrixstr
	Note/k cmp_mass_Neg "This wave contains the negative values in wave cmp_mass"
	Note/k cmp_mass_Neg_indx "This wave contains the index number of negative values in wave cmp_mass"
End

//This function finds the same index number in two waves and extract them out
Function CrossComparingIndx()

	string wave1str="TooNeg_indx"
	string wave2str="cmp_mass_Neg_indx"
	string newwavestr = "Neg2Remove_indx"
	
	prompt wave1str, "Select first index wave", popup WaveList("*",";","")
	prompt wave2str, "Select second index wave", popup WaveList("*",";","")
	prompt newwavestr, "Name of the newly created index wave"
	Doprompt "Select waves", wave1str,wave2str, newwavestr
	
	if(V_flag)
		abort
	endif
	
	wave wave1=$wave1str
	wave wave2=$wave2str
	
	variable npnts1=numpnts(wave1)
	variable npnts2=numpnts(wave2)
	
	variable i,j,count
	count=0
	
	make/o/d/n=0 $(newwavestr)
	wave newwave=$newwavestr
	
	for(i=0;i<npnts1;i+=1)
		for(j=0;j<npnts2;j+=1)
			if(wave2[j]==wave1[i])
				InsertPoints count,1, newwave
				newwave[count] = wave1[i]
				count+=1
			endif
		endfor
	endfor
	
	Note/k newwave "This wave contains the values that are found both in wave "+nameofwave(wave1)+" and wave "+nameofwave(wave2)
	
End

//This function combines two index waves to include the index number in both waves
Function CombineIndx()

	string wave1str="TooNeg_indx"
	string wave2str="cmp_mass_Neg_indx"
	string newwavestr = "Neg2Remove_indx"
	
	prompt wave1str, "Select first index wave", popup WaveList("*",";","")
	prompt wave2str, "Select second index wave", popup WaveList("*",";","")
	prompt newwavestr, "Name of the newly created index wave"
	Doprompt "Select waves", wave1str,wave2str, newwavestr
	
	if(V_flag)
		abort
	endif
	
	wave wave1=$wave1str
	wave wave2=$wave2str
	
	variable npnts1=numpnts(wave1)
	variable npnts2=numpnts(wave2)
	
	variable i,j,count
	count=0
	
	Concatenate/NP/O {wave1,wave2},$newwavestr
	wave newwave=$newwavestr

	Sort newwave,newwave
	
	variable npnts=numpnts(newwave)
	variable tempVal = newwave[0]
	
	for(i=1;i<npnts;i+=1)
		if(newwave[i]==tempVal)
			newwave[i] = 10000
		else
			tempVal= newwave[i]
		endif		
	endfor
	
	i=0
	
	do
		if(newwave[i]==10000)
			Deletepoints i,1,newwave
		endif
		i+=1
	while(i<numpnts(newwave))
	
	Note/k newwave "This wave contains the values combined from wave "+nameofwave(wave1)+" and wave "+nameofwave(wave2)
	
End

Function CombineIndx_np(wave1,wave2,newwavestr)

	wave wave1, wave2
	string newwavestr
	
	variable npnts1=numpnts(wave1)
	variable npnts2=numpnts(wave2)
	
	variable i,j,count
	count=0
	
	Concatenate/NP/O {wave1,wave2},$newwavestr
	wave newwave=$newwavestr

	Sort newwave,newwave
	
	variable npnts=numpnts(newwave)
	variable tempVal = newwave[0]
	
	for(i=1;i<npnts;i+=1)
		if(newwave[i]==tempVal)
			newwave[i] = 10000
		else
			tempVal= newwave[i]
		endif		
	endfor
	
	i=0
	
	do
		if(newwave[i]==10000)
			Deletepoints i,1,newwave
		endif
		i+=1
	while(i<numpnts(newwave))
	
	Note/k newwave "This wave contains the values combined from wave "+nameofwave(wave1)+" and wave "+nameofwave(wave2)
	
End



//This function aims to remove the points/columns from 1D/2D wave according to the index in the indexwave
Function RemoveWeirdIons()
	string wavestr="cmp_matrix"
	string indxwavestr="Neg2remove_indx"
	
	prompt wavestr, "Select wave to remove points", popup WaveList("*",";","")
	prompt indxwavestr, "Select index wave", popup WaveList("*",";","")
	Doprompt "Select waves",wavestr,indxwavestr
	
	if(V_flag)
		abort
	endif

	wave indxwave=$indxwavestr	
	variable nindx=numpnts(indxwave)
	
	variable i
	
	if(Wavetype($wavestr,1)==1)
	
		wave wave2remove = $wavestr
		variable nrows=dimsize(wave2remove,0)
		variable ncols=dimsize(wave2remove,1)
		make/o/d/n=(nrows,ncols) $(wavestr+"_rf") //_rf denotes for refined
		wave refinedwave= $(wavestr+"_rf")
		refinedwave = wave2remove
		Note/k refinedwave "This is the refined wave of "+wavestr+" and weird points/columns are removed according to index in wave "+indxwavestr

		if(ncols==0)
		//1D wave
			for(i=nindx-1;i>=0;i-=1)
			//for(i=0;i<nindx;i+=1) <- this is wrong. must delete backwards
				DeletePoints indxwave[i],1,refinedwave
			endfor
		else
		//2D wave
			for(i=nindx-1;i>=0;i-=1)
			//for(i=0;i<nindx;i+=1) <- this is wrong. must delete backwards
				DeletePoints/M=1 indxwave[i],1,refinedwave
			endfor
		endif
	
	elseif(Wavetype($wavestr,1)==2)
	
		wave/t wave2removet = $wavestr
		variable nrowst=dimsize(wave2removet,0)
		make/o/t/n=(nrowst) $(wavestr+"_rf")
		wave/t refinedwavet= $(wavestr+"_rf")
		for(i=0;i<nrowst;i+=1)
			refinedwavet[i] = wave2removet[i]
		endfor
		Note/k refinedwavet "This is the refined wave of "+wavestr+" and weird points/columns are removed according to index in wave "+indxwavestr

		for(i=nindx-1;i>=0;i-=1)
			DeletePoints indxwave[i],1,refinedwavet
		endfor
		
	endif
	
	
	
	
End


Function RemoveWeirdIons_np(ywave,indxwave)

	wave ywave, indxwave

	variable nindx=numpnts(indxwave)
	
	string wavestr = nameofwave(ywave)
	string indxwavestr = nameofwave(indxwave)
	
	variable i
	
	if(Wavetype($wavestr,1)==1)
	
		wave wave2remove = $wavestr
		variable nrows=dimsize(wave2remove,0)
		variable ncols=dimsize(wave2remove,1)
		make/o/d/n=(nrows,ncols) $(wavestr+"_rf") //_rf denotes for refined
		wave refinedwave= $(wavestr+"_rf")
		refinedwave = wave2remove
		Note/k refinedwave "This is the refined wave of "+wavestr+" and weird points/columns are removed according to index in wave "+indxwavestr

		if(ncols==0)
		//1D wave
			for(i=nindx-1;i>=0;i-=1)
			//for(i=0;i<nindx;i+=1) <- this is wrong. must delete backwards
				DeletePoints indxwave[i],1,refinedwave
			endfor
		else
		//2D wave
			for(i=nindx-1;i>=0;i-=1)
			//for(i=0;i<nindx;i+=1) <- this is wrong. must delete backwards
				DeletePoints/M=1 indxwave[i],1,refinedwave
			endfor
		endif
	
	elseif(Wavetype($wavestr,1)==2)
	
		wave/t wave2removet = $wavestr
		variable nrowst=dimsize(wave2removet,0)
		make/o/t/n=(nrowst) $(wavestr+"_rf")
		wave/t refinedwavet= $(wavestr+"_rf")
		for(i=0;i<nrowst;i+=1)
			refinedwavet[i] = wave2removet[i]
		endfor
		Note/k refinedwavet "This is the refined wave of "+wavestr+" and weird points/columns are removed according to index in wave "+indxwavestr

		for(i=nindx-1;i>=0;i-=1)
			DeletePoints indxwave[i],1,refinedwavet
		endfor
		
	endif
	
	
	
	
End



//This function graphs columns in cmp_matrix specified by indexwave as different panels in one graph
Function GraphIndx_col()
	string ywvstr="cmp_matrix"
	string xwvstr="time_s"
	string indxwavestr="TooNeg_indx"
	variable ncols=3
	
	prompt ywvstr, "Select matrix for y axis", popup WaveList("*",";","")
	prompt xwvstr, "Select wave for x axis", popup WaveList("*",";","")
	prompt indxwavestr, "Select index wave", popup WaveList("*",";","")
	prompt ncols, "Columns of panels of graph"
	Doprompt "Select waves",ywvstr,xwvstr,indxwavestr,ncols
	
	if(V_flag)
		abort
	endif
	
	wave ywv=$ywvstr
	wave xwv=$xwvstr
	wave indxwave=$indxwavestr
	
	variable nwvs = numpnts(indxwave) // number of waves
	variable nrows = ceil(nwvs/ncols) // number of rows
	
	variable col_offset = 0.05
	variable col_len 
	variable col_space 
	variable col_multiplier = 2
	
	col_offset = 0.05
	//Let col_len=col_multiplier*col_space, then col_offset+ncols*col_len+(ncols-1)*col_space=1
	//so that
	col_space = (1-col_offset)/((1+col_multiplier)*ncols-1)
	col_len=col_multiplier*col_space
	
	variable row_len
	//if(nrows>1)
		//row_len = (1-0.1/(nrows/2))*1/nrows
	//else
		//row_len = 1
	//endif
	//variable row_space = (1-row_len*nrows)/3
	variable row_space = 0.1/nrows
	row_len = 1/nrows
	
	//print row_len,row_space,col_len,col_space
	
	string yaxstr, xaxstr
	variable i, j
	variable wvcount = 0
	variable index
	
	variable xlowFra,xhighFra,ylowFra,yhighFra
	string abcstr="a;b;c;d;e;f;g;h;i;j;k;l;m;n;o;p;q;r;s;t;u;v;w;x;y;z"
	string NameofPanel
	
	string graphname = "TG_"+indxwavestr
	display/N=$graphname as graphname
	modifygraph width=500, height=500
	for(i=0;i<nrows;i+=1)
		for(j=0;j<ncols;j+=1)
			if(wvcount > nwvs-1)
				abort
			endif
			index=j+1+ncols*i
			yaxstr = "y" + num2istr(index)
			xaxstr = "x" + num2istr(index)
			appendtograph/l=$yaxstr/b=$xaxstr ywv[][indxwave[wvcount]] vs xwv
//			print 1-row_len*(i+1)-row_space*i,1-(row_len)*i
//			ModifyGraph axisEnab($yaxstr)={1-row_len*(i+1),1-(row_len+row_space)*i},axisEnab($xaxstr)={col_len*j+col_space*j,(col_len)*(j+1)}
			// still need to work on x spacing
//			ModifyGraph axisEnab($yaxstr)={1-row_len*(i+1)+row_space*(nrows-i-1),1-(row_len)*i-row_space*(i)},axisEnab($xaxstr)={col_len*j+col_space*j,col_len2*(j+1)-col_space*(ncols-j-1)}
			if(j==0)
				//ModifyGraph axisEnab($yaxstr)={1-row_len*(i+1)+row_space*(nrows-i),1-(row_len)*i-row_space*(i)},axisEnab($xaxstr)={col_offset+col_len*j,col_offset+col_len*(j+1)} // col_len*j+col_space*j,col_len2*(j+1)-col_space*(ncols-j-1)}
				ylowFra = 1-row_len*(i+1)+row_space //*(nrows-i)
				yhighFra = 1-(row_len)*i-row_space //*(i)
				xlowFra = col_offset+col_len*j
				xhighFra = col_offset+col_len*(j+1)
				ModifyGraph axisEnab($yaxstr)={ylowFra,yhighFra},axisEnab($xaxstr)={xlowFra,xhighFra}
			else
				ylowFra = 1-row_len*(i+1)+row_space //*(nrows-i)
				yhighFra = 1-(row_len)*i-row_space //*(i)
				xlowFra = col_offset+col_len*j+col_space*j
				xhighFra = col_offset+col_len*(j+1)+col_space*j
				ModifyGraph axisEnab($yaxstr)={ylowFra,yhighFra},axisEnab($xaxstr)={xlowFra,xhighFra}
			endif
			//print ylowFra,yhighFra,xlowFra,xhighFra			
			ModifyGraph freePos($yaxstr)={0,$xaxstr},freePos($xaxstr)={0,$yaxstr}
			ModifyGraph noLabel($xaxstr)=1//, noLabel($yaxstr)=1
			//NameofPanel= "\\Z16"+stringfromlist(index-1,abcstr,";")+")"
			NameofPanel=num2istr(indxwave[wvcount])
			TextBox/C/N=$("name"+num2istr(index))/A=MC/X=(xhighFra*100-50)/Y=(yhighFra*100-50) nameofpanel
			wvcount +=1
		endfor
	endfor
	
End

Function GraphIndx_col_np(ywv,xwv,indxwave,ncols)
	
	wave ywv
	wave xwv
	wave indxwave
	variable ncols
	
	string ywvstr=nameofwave(ywv)
	string xwvstr=nameofwave(xwv)
	string indxwavestr=nameofwave(indxwave)
	
	variable nwvs = numpnts(indxwave) // number of waves
	variable nrows = ceil(nwvs/ncols) // number of rows
	
	variable col_offset
	variable col_len 
	variable col_space 
	variable col_multiplier = 2
	
	col_offset = 0.05
	//Let col_len=col_multiplier*col_space, then col_offset+ncols*col_len+(ncols-1)*col_space=1
	//so that
	col_space = (1-col_offset)/((1+col_multiplier)*ncols-1)
	col_len=col_multiplier*col_space
	
	variable row_len
	//if(nrows>1)
		//row_len = (1-0.1/(nrows/2))*1/nrows
	//else
		//row_len = 1
	//endif
	//variable row_space = (1-row_len*nrows)/3
	variable row_space = 0.1/nrows
	row_len = 1/nrows
	
	//print row_len,row_space,col_len,col_space
	
	string yaxstr, xaxstr
	variable i, j
	variable wvcount = 0
	variable index
	
	variable xlowFra,xhighFra,ylowFra,yhighFra
	string abcstr="a;b;c;d;e;f;g;h;i;j;k;l;m;n;o;p;q;r;s;t;u;v;w;x;y;z"
	string NameofPanel
	
	string graphname =  "TG_"+indxwavestr
	//display
	display/N=$graphname as graphname
	modifygraph width=500, height=500
	
	for(i=0;i<nrows;i+=1)
		for(j=0;j<ncols;j+=1)
			if(wvcount > nwvs-1)
				return 1
			endif
			index=j+1+ncols*i
			yaxstr = "y" + num2istr(index)
			xaxstr = "x" + num2istr(index)
			appendtograph/l=$yaxstr/b=$xaxstr ywv[][indxwave[wvcount]] vs xwv
//			print 1-row_len*(i+1)-row_space*i,1-(row_len)*i
//			ModifyGraph axisEnab($yaxstr)={1-row_len*(i+1),1-(row_len+row_space)*i},axisEnab($xaxstr)={col_len*j+col_space*j,(col_len)*(j+1)}
			// still need to work on x spacing
//			ModifyGraph axisEnab($yaxstr)={1-row_len*(i+1)+row_space*(nrows-i-1),1-(row_len)*i-row_space*(i)},axisEnab($xaxstr)={col_len*j+col_space*j,col_len2*(j+1)-col_space*(ncols-j-1)}
			if(j==0)
				//ModifyGraph axisEnab($yaxstr)={1-row_len*(i+1)+row_space*(nrows-i),1-(row_len)*i-row_space*(i)},axisEnab($xaxstr)={col_offset+col_len*j,col_offset+col_len*(j+1)} // col_len*j+col_space*j,col_len2*(j+1)-col_space*(ncols-j-1)}
				ylowFra = 1-row_len*(i+1)+row_space //*(nrows-i)
				yhighFra = 1-(row_len)*i-row_space //*(i)
				xlowFra = col_offset+col_len*j
				xhighFra = col_offset+col_len*(j+1)
				ModifyGraph axisEnab($yaxstr)={ylowFra,yhighFra},axisEnab($xaxstr)={xlowFra,xhighFra}
			else
				ylowFra = 1-row_len*(i+1)+row_space //*(nrows-i)
				yhighFra = 1-(row_len)*i-row_space //*(i)
				xlowFra = col_offset+col_len*j+col_space*j
				xhighFra = col_offset+col_len*(j+1)+col_space*j
				ModifyGraph axisEnab($yaxstr)={ylowFra,yhighFra},axisEnab($xaxstr)={xlowFra,xhighFra}
			endif
			//print ylowFra,yhighFra,xlowFra,xhighFra			
			ModifyGraph freePos($yaxstr)={0,$xaxstr},freePos($xaxstr)={0,$yaxstr}
			ModifyGraph noLabel($xaxstr)=1//, noLabel($yaxstr)=1
			//NameofPanel= "\\Z16"+stringfromlist(index-1,abcstr,";")+")"
			NameofPanel=num2istr(indxwave[wvcount])
			TextBox/C/N=$("name"+num2istr(index))/A=MC/X=(xhighFra*100-50)/Y=(yhighFra*100-50) nameofpanel
			wvcount +=1
		endfor
	endfor

End



//This function plots the mass spectrum and highlights the specific ions according to index wave
Function GraphIndx_massSpec()

	string ywvstr="cmp_mass"
	string indxwavestr="Neg2remove_indx"
	
	prompt ywvstr, "Select wave of mass spectrum", popup WaveList("*",";","")
	prompt indxwavestr, "Select index wave", popup WaveList("*",";","")
	Doprompt "Select waves",ywvstr,indxwavestr
	
	if(V_flag)
		abort
	endif
	
	wave ywv=$(ywvstr)
	wave indxwave=$(indxwavestr)
	
	string graphname = "MassSpectrum_"+indxwavestr
	
	Display/N=$graphname as graphname
	
	appendtograph ywv
	
	variable nrows = numpnts(ywv)
	variable nindx = numpnts(indxwave)
	
	string ywvweirdstr = ywvstr+"_weird"
	make/o/d/n=(nrows) $(ywvweirdstr)
	wave ywvweird = $(ywvweirdstr)
	
	Note/k ywvweird "This wave is nan except for the points indicated by index wave "+indxwavestr+" have the same values as in wave "+ywvstr
	
	ywvweird = nan
	
	variable i
	for(i=0;i<nindx;i+=1)
		ywvweird[indxwave[i]] = ywv[indxwave[i]]
	endfor
	
	appendtograph ywvweird
	
	ModifyGraph mode($ywvstr)=1,rgb($ywvstr)=(65535,32768,32768)
	ModifyGraph mode($ywvweirdstr)=8,marker($ywvweirdstr)=19,rgb($ywvweirdstr)=(0,0,0)
	ModifyGraph lsize=1.5
	
	ModifyGraph fSize=12
	Label left "Absolute mass"
	Label bottom "Compound index"

	modifygraph width=500, height=120
	
End

Function GraphIndx_massSpec_np(ywv,indxwave)
	wave ywv, indxwave

	string ywvstr=nameofwave(ywv)
	string indxwavestr=nameofwave(indxwave)
	
	string graphname = "MassSpectrum_"+indxwavestr
	
	Display/N=$graphname as graphname
	//display
	appendtograph ywv
	
	variable nrows = numpnts(ywv)
	variable nindx = numpnts(indxwave)
	
	string ywvweirdstr = ywvstr+"_weird"
	make/o/d/n=(nrows) $(ywvweirdstr)
	wave ywvweird = $(ywvweirdstr)
	
	Note/k ywvweird "This wave is nan except for the points indicated by index wave "+indxwavestr+" have the same values as in wave "+ywvstr
	
	ywvweird = nan
	
	variable i
	for(i=0;i<nindx;i+=1)
		ywvweird[indxwave[i]] = ywv[indxwave[i]]
	endfor
	
	appendtograph ywvweird
	
	ModifyGraph mode($ywvstr)=1,rgb($ywvstr)=(65535,32768,32768)
	ModifyGraph mode($ywvweirdstr)=8,marker($ywvweirdstr)=19,rgb($ywvweirdstr)=(0,0,0)
	ModifyGraph lsize=1.5
	
	ModifyGraph fSize=12
	Label left "Absolute mass"
	Label bottom "Compound index"

	modifygraph width=500, height=120
	
End


// This function smoothes the thermogram of every compound by using moving average of X points
//Thermograms are represented in columns of the matrix
Function MovingAverageMatrix_Col()
	
	string ymatrixstr = "cmp_matrix"
	variable avgpnt = 35
	
	prompt ymatrixstr, "Select matrix to smooth", popup WaveList("*",";","")
	prompt avgpnt, "Number of points used for moving average"
	Doprompt "Select waves",ymatrixstr, avgpnt
	
	if(V_flag)
		abort
	endif
	
	wave ymatrix=$ymatrixstr
	
	variable nrows=dimsize(ymatrix,0)
	variable ncols=dimsize(ymatrix,1)
	
	string newymatrix_str=nameofwave(ymatrix)+"_Smt"+num2str(avgpnt)
	
	make/o/d/n=(nrows)/Free cmp_thermo
	make/o/d/n=(nrows-avgpnt+1,ncols) $newymatrix_str
	make/o/d/n=(avgpnt)/Free tempwave
	wave newymatrix=$newymatrix_str
	
	variable i,j,k
	
	for(i=0;i<ncols;i+=1)
		cmp_thermo=ymatrix[p][i]
		for(j=avgpnt-1;j<nrows;j+=1)
			for(k=0;k<avgpnt;k+=1)
				tempwave[k]=cmp_thermo[j-avgpnt+1+k]
			endfor
			wavestats/q tempwave
			newymatrix[j-avgpnt+1][i]=V_avg
		endfor
	endfor
	
	Note/k newymatrix "This matrix stores the smoothed curve of every column of the matrix ' "+nameofwave(ymatrix)+" ' "
	Note newymatrix "The smoothing is done by using "+num2str(avgpnt)+" past points moving average"
	
End

Function MovingAverageMatrix_Col_np(ymatrix,avgpnt)
	
	wave ymatrix
	variable avgpnt
	
	string ymatrixstr=nameofwave(ymatrix)
	
	variable nrows=dimsize(ymatrix,0)
	variable ncols=dimsize(ymatrix,1)
	
	string newymatrix_str=nameofwave(ymatrix)+"_Smt"+num2str(avgpnt)
	
	make/o/d/n=(nrows)/Free cmp_thermo
	make/o/d/n=(nrows-avgpnt+1,ncols) $newymatrix_str
	make/o/d/n=(avgpnt)/Free tempwave
	wave newymatrix=$newymatrix_str
	
	variable i,j,k
	
	for(i=0;i<ncols;i+=1)
		cmp_thermo=ymatrix[p][i]
		for(j=avgpnt-1;j<nrows;j+=1)
			for(k=0;k<avgpnt;k+=1)
				tempwave[k]=cmp_thermo[j-avgpnt+1+k]
			endfor
			wavestats/q tempwave
			newymatrix[j-avgpnt+1][i]=V_avg
		endfor
	endfor
	
	Note/k newymatrix "This matrix stores the smoothed curve of every column of the matrix ' "+nameofwave(ymatrix)+" ' "
	Note newymatrix "The smoothing is done by using "+num2str(avgpnt)+" past points moving average"
	
End


////////////////////////////////////////////////////////////////////
Function MovingAverageWave()
	
	string xwavestr = "time_s"
	variable avgpnt = 35
	
	prompt xwavestr, "Select wave to smooth", popup WaveList("*",";","")
	prompt avgpnt, "Number of points used for moving average"
	Doprompt "Select waves",xwavestr, avgpnt
	
	if(V_flag)
		abort
	endif
	
	wave xwave = $xwavestr
	
	variable nrows=numpnts(xwave)
	
	string newxwave_str=nameofwave(xwave)+"_Smt"+num2str(avgpnt)
	
	make/o/d/n=(nrows-avgpnt+1) $newxwave_Str
	wave newxwave=$newxwave_str
	
	make/o/d/n=(avgpnt)/Free tempwave
	variable i,k
	for(i=avgpnt-1;i<nrows;i+=1)
		for(k=0;k<avgpnt;k+=1)
			tempwave[k]=xwave[i-avgpnt+1+k]
		endfor
		wavestats/q tempwave
		newxwave[i-avgpnt+1]=V_avg
	endfor
	
	Note/k newxwave "This wave stores the "+num2str(avgpnt)+" points moving average of wave "+nameofwave(xwave)

	
End

Function MovingAverageWave_np(xwave, avgpnt)
	
	wave xwave
	variable avgpnt
	
	string xwavestr=nameofwave(xwave)
	
	variable nrows=numpnts(xwave)
	
	string newxwave_str=nameofwave(xwave)+"_Smt"+num2str(avgpnt)
	
	make/o/d/n=(nrows-avgpnt+1) $newxwave_Str
	wave newxwave=$newxwave_str
	
	make/o/d/n=(avgpnt)/Free tempwave
	variable i,k
	for(i=avgpnt-1;i<nrows;i+=1)
		for(k=0;k<avgpnt;k+=1)
			tempwave[k]=xwave[i-avgpnt+1+k]
		endfor
		wavestats/q tempwave
		newxwave[i-avgpnt+1]=V_avg
	endfor
	
	Note/k newxwave "This wave stores the "+num2str(avgpnt)+" points moving average of wave "+nameofwave(xwave)

	
End



//Normalize the columns in the matrix to the smoothed maximum in each column
Function Normalize2SmoothedMax()
	
	string smtymatrixstr = "cmp_matrix_rf_smt35"
	string ymatrixstr = "cmp_matrix_rf"
	
	prompt smtymatrixstr, "Select smoothed matrix", popup WaveList("*",";","")
	prompt ymatrixstr, "Select unsmoothed matrix", popup WaveList("*",";","")
	Doprompt "Select waves",smtymatrixstr,ymatrixstr
	
	if(V_flag)
		abort
	endif
	
	wave smtymatrix = $(smtymatrixstr)
	wave ymatrix = $(ymatrixstr)
	
	variable nrows=dimsize(ymatrix,0)
	variable ncols=dimsize(ymatrix,1)
	variable nrowssmt=dimsize(smtymatrix,0)
	variable ncolssmt=dimsize(smtymatrix,1)
	
	if(ncols!=ncolssmt)
		print "please make sure two matrices have the same number of columns"
		abort
	endif
	
	make/o/d/n=(nrows,ncols) $(nameofwave(ymatrix)+"_smtnorm")
	wave newymatrixsmt=$(nameofwave(ymatrix)+"_smtnorm")
	make/o/d/n=(nrows,ncols) $(nameofwave(ymatrix)+"_norm")
	wave newymatrix=$(nameofwave(ymatrix)+"_norm")

	make/o/d/n=(nrowssmt)/Free tempwave
	make/o/d/n=(nrows)/Free tempwave1
	
	variable i
	variable SmtMaxT
	variable SmtMaxIndx
	
	for(i=0;i<ncols;i+=1)
		tempwave=smtymatrix[p][i]
		wavestats/q tempwave
		tempwave1=ymatrix[p][i]/V_max
		newymatrixsmt[][i]=tempwave1[p]
		tempwave=ymatrix[p][i]
		wavestats/q tempwave
		tempwave1=ymatrix[p][i]/V_max
		newymatrix[][i]=tempwave1[p]
	endfor
	
	Note/k newymatrixsmt "Every column in matrix "+ymatrixstr+" is normalized to the maximum in the corresponding column in matrix "+smtymatrixstr
	Note/k newymatrix "Every column in matrix "+ymatrixstr+" is normalized to the maximum in the column"
	
End

Function Normalize2SmoothedMax_np(smtymatrix,ymatrix)
	
	wave smtymatrix
	wave ymatrix
	
	string smtymatrixstr = nameofwave(smtymatrix)
	string ymatrixstr = nameofwave(ymatrix)
	
	variable nrows=dimsize(ymatrix,0)
	variable ncols=dimsize(ymatrix,1)
	variable nrowssmt=dimsize(smtymatrix,0)
	variable ncolssmt=dimsize(smtymatrix,1)
	
	if(ncols!=ncolssmt)
		print "please make sure two matrices have the same number of columns"
		abort
	endif
	
	make/o/d/n=(nrows,ncols) $(nameofwave(ymatrix)+"_smtnorm")
	wave newymatrixsmt=$(nameofwave(ymatrix)+"_smtnorm")
	make/o/d/n=(nrows,ncols) $(nameofwave(ymatrix)+"_norm")
	wave newymatrix=$(nameofwave(ymatrix)+"_norm")

	make/o/d/n=(nrowssmt)/Free tempwave
	make/o/d/n=(nrows)/Free tempwave1
	
	variable i
	variable SmtMaxT
	variable SmtMaxIndx
	
	for(i=0;i<ncols;i+=1)
		tempwave=smtymatrix[p][i]
		wavestats/q tempwave
		tempwave1=ymatrix[p][i]/V_max
		newymatrixsmt[][i]=tempwave1[p]
		tempwave=ymatrix[p][i]
		wavestats/q tempwave
		tempwave1=ymatrix[p][i]/V_max
		newymatrix[][i]=tempwave1[p]
	endfor
	
	Note/k newymatrixsmt "Every column in matrix "+ymatrixstr+" is normalized to the maximum in the corresponding column in matrix "+smtymatrixstr
	Note/k newymatrix "Every column in matrix "+ymatrixstr+" is normalized to the maximum in the column"
	
End


//This function averages every x points for the range from startpnt to endpnt
//This function serves for down-weighting certain part of the data
Function AverageWave_Subrange()
	
	string ywavestr = "time_s"
	variable avgpnt = 10
	variable startpnt = 240
	variable endpnt = 800
	string suffix = "_10a"
	
	prompt ywavestr, "Select wave to average partly", popup WaveList("*",";","")
	prompt avgpnt, "Number of points to average "
	prompt startpnt, "Start pnt of the subrange"
	prompt endpnt, "End pnt of the subrange"
	prompt suffix, "Suffix to add to the new wave"
	Doprompt "Select waves",ywavestr, avgpnt, startpnt, endpnt, suffix
	
	if(V_flag)
		abort
	endif

	wave ywave = $ywavestr
	
	variable nrows=numpnts(ywave)
	variable nrowspart=ceil((endpnt-startpnt+1)/avgpnt)
	variable nrowsnew=nrows-(endpnt-startpnt+1)+nrowspart
	
	make/o/d/n=(avgpnt)/Free tempwave
	make/o/d/n=(nrowsnew) $(ywavestr+suffix)
	wave newywave=$(ywavestr+suffix)

	variable i,j,newi
	newi=startpnt
	
	for(i=0;i<startpnt;i+=1)
		newywave[i]=ywave[i]
	endfor
	
	for(i=startpnt;i<=endpnt;i+=avgpnt)
		tempwave=nan
		if(i+avgpnt>endpnt)
			for(j=0;j<endpnt-i+1;j+=1)
				tempwave[j]=ywave[i+j]
			endfor
			wavestats/q tempwave
			newywave[newi]=V_avg
		else
			for(j=0;j<avgpnt;j+=1)
				tempwave[j]=ywave[i+j]
			endfor
			wavestats/q tempwave
			newywave[newi]=V_avg
		endif
		newi+=1
	endfor
	
	for(i=endpnt+1;i<nrows;i+=1)
		newywave[newi]=ywave[i]
		newi+=1
	endfor
	
	string notestr="This wave stores the data in wave "+ywavestr+", but data from point "+num2str(startpnt)+" and "+num2str(endpnt)+" are averaged by every "+num2str(avgpnt)+" points" 
	Note/k newywave notestr
	
End


Function AverageWave_Subrange_np(ywave,avgpnt,startpnt,endpnt,suffix)
	wave ywave
	variable avgpnt, startpnt, endpnt
	string suffix
	
	variable nrows=numpnts(ywave)
	variable nrowspart=ceil((endpnt-startpnt+1)/avgpnt)
	variable nrowsnew=nrows-(endpnt-startpnt+1)+nrowspart
	
	string ywavestr=nameofwave(ywave)
	
	make/o/d/n=(avgpnt)/Free tempwave
	make/o/d/n=(nrowsnew) $(ywavestr+suffix)
	wave newywave=$(ywavestr+suffix)

	variable i,j,newi
	newi=startpnt
	
	for(i=0;i<startpnt;i+=1)
		newywave[i]=ywave[i]
	endfor
	
	for(i=startpnt;i<=endpnt;i+=avgpnt)
		tempwave=nan
		if(i+avgpnt>endpnt)
			for(j=0;j<endpnt-i+1;j+=1)
				tempwave[j]=ywave[i+j]
			endfor
			wavestats/q tempwave
			newywave[newi]=V_avg
		else
			for(j=0;j<avgpnt;j+=1)
				tempwave[j]=ywave[i+j]
			endfor
			wavestats/q tempwave
			newywave[newi]=V_avg
		endif
		newi+=1
	endfor
	
	for(i=endpnt+1;i<nrows;i+=1)
		newywave[newi]=ywave[i]
		newi+=1
	endfor
	
	string notestr="This wave stores the data in wave "+ywavestr+", but data from point "+num2str(startpnt)+" and "+num2str(endpnt)+" are averaged by every "+num2str(avgpnt)+" points" 
	Note/k newywave notestr
	
End


//This function assumes that the average is done among rows, so that the number of columns of the new matrix stays the same
//This function serves for down-weighting certain part of the data
Function AverageMatrix_Subrange()
	
	string ywavestr = "cmp_matrix_rf_smtnorm"
	variable avgpnt = 10
	variable startpnt = 241
	variable endpnt = 831
	string suffix = "_10a"
	
	prompt ywavestr, "Select matrix to average partly", popup WaveList("*",";","")
	prompt avgpnt, "Number of points to average"
	prompt startpnt, "Start pnt of the subrange"
	prompt endpnt, "End pnt of the subrange"
	prompt suffix, "Suffix to add to the new waves"
	Doprompt "Select waves",ywavestr, avgpnt, startpnt, endpnt, suffix
	
	if(V_flag)
		abort
	endif
	
	wave ywave = $ywavestr
	
	variable nrows=dimsize(ywave,0)
	variable ncols=dimsize(ywave,1)
	variable nrowspart=ceil((endpnt-startpnt+1)/avgpnt)
	variable nrowsnew=nrows-(endpnt-startpnt+1)+nrowspart
	
	make/o/d/n=(nrows) tempywave
	make/o/d/n=(nrowsnew,ncols) $(nameofwave(ywave)+suffix)
	wave newywave=$(nameofwave(ywave)+suffix)
	
	variable i,j,newi
	newi=startpnt
	
	for(i=0;i<ncols;i+=1)
		tempywave=ywave[p][i]
		Averagewave_Subrange_np(tempywave,avgpnt, startpnt,endpnt,suffix)
		wave tempnewywave=$("tempywave"+suffix)
		newywave[][i]=tempnewywave[p]
	endfor
	
	killwaves tempywave,$("tempywave"+suffix)
	
	string notestr="This matrix stores the data in matrix "+ywavestr+", but data from row "+num2str(startpnt)+" and "+num2str(endpnt)+" are averaged by every "+num2str(avgpnt)+" rows" 
	Note/k newywave notestr
	
End


Function AverageMatrix_Subrange_np(ywave,avgpnt,startpnt,endpnt,suffix)
	wave ywave
	variable avgpnt, startpnt, endpnt
	string suffix
	
	string ywavestr = nameofwave(ywave)
	
	variable nrows=dimsize(ywave,0)
	variable ncols=dimsize(ywave,1)
	variable nrowspart=ceil((endpnt-startpnt+1)/avgpnt)
	variable nrowsnew=nrows-(endpnt-startpnt+1)+nrowspart
	
	make/o/d/n=(nrows) tempywave
	make/o/d/n=(nrowsnew,ncols) $(nameofwave(ywave)+suffix)
	wave newywave=$(nameofwave(ywave)+suffix)
	
	variable i,j,newi
	newi=startpnt
	
	for(i=0;i<ncols;i+=1)
		tempywave=ywave[p][i]
		Averagewave_Subrange_np(tempywave,avgpnt, startpnt,endpnt,suffix)
		wave tempnewywave=$("tempywave"+suffix)
		newywave[][i]=tempnewywave[p]
	endfor
	
	killwaves tempywave,$("tempywave"+suffix)
	
	string notestr="This matrix stores the data in matrix "+ywavestr+", but data from row "+num2str(startpnt)+" and "+num2str(endpnt)+" are averaged by every "+num2str(avgpnt)+" rows" 
	Note/k newywave notestr
	
End

//This function calculates the difference between unsmoothed and smoothed curve (column in a matrix)
//The subrange is based on the smoothed waves. Note that the smoothing is done by moving average, so the number of rows for smoothed and unsmoothed thermograms are different
Function CalculateNoise_Subrange()

	string ymatrixstr = "cmp_matrix_rf_smtnorm"
	string ymatrixsmtstr = "cmp_matrix_rf_smtnorm_smt35"
	string xwavestr = "temp_c"
	string xwavesmtstr = "temp_c_smt35"
	variable startpnt = 0
	variable endpnt = 223
	string suffix = "_ar35"
	
	prompt ymatrixstr, "Select unsmoothed matrix", popup WaveList("*",";","")
	prompt ymatrixsmtstr, "Select smoothed matrix", popup WaveList("*",";","")
	prompt xwavestr, "Select unsmoothed T wave", popup WaveList("*",";","")
	prompt xwavesmtstr, "Select smoothed T wave", popup WaveList("*",";","")
	prompt startpnt, "Start pnt of the subrange"
	prompt endpnt, "End pnt of the subrange"
	prompt suffix, "Suffix to add to the new wave"
	Doprompt "Select waves",ymatrixstr, ymatrixsmtstr,xwavestr,xwavesmtstr, startpnt, endpnt, suffix
	
	if(V_flag)
		abort
	endif
	
	wave ymatrix = $ymatrixstr
	wave smt_ymatrix = $ymatrixsmtstr
	wave xwave = $xwavestr
	wave smt_xwave = $xwavesmtstr
	
	variable nrows=dimsize(ymatrix,0)
	variable ncols=dimsize(ymatrix,1)
	variable nrows_smt=dimsize(smt_ymatrix,0)
	variable ncols_smt=dimsize(smt_ymatrix,1)
	
	if(ncols!=ncols_smt)
		print "Please make sure two matrices have the same number of columns"
		abort
	endif
	
	make/o/d/n=(nrows_smt,ncols_smt)/Free newymatrix
	make/o/d/n=(ncols_smt) $("AvgNoise"+suffix)
	make/o/d/n=(ncols_smt) $("NoiseSDV"+suffix)
	make/o/d/n=(nrows_smt)/Free tempwave
	
	wave newwave=$("AvgNoise"+suffix)
	wave newwave_sdv=$("NoiseSDV"+suffix)
	
	variable i,j
	
	for(i=0;i<ncols;i+=1)
		for(j=startpnt;j<endpnt;j+=1)
			findlevel/Q/P xwave,smt_xwave[j]
			newymatrix[j][i]=abs(ymatrix(V_levelx)[i]-smt_ymatrix[j][i])
		endfor
		tempwave=newymatrix[p][i]
		wavestats/q tempwave
		newwave[i]=V_avg
		newwave_sdv[i]=V_sdev
		
	endfor
	
	Note/k newwave "This wave stores the average difference between smoothed ( from matrix "+ymatrixstr+" ) and unsmoothed (from "+ymatrixsmtstr+" ) as noise level from row "+num2str(startpnt)+" to row "+num2str(endpnt)
	Note/k newwave_sdv "This wave stores the standard deviation of noise. Calculated along with the wave "+nameofwave(newwave)
End



Function CalculateNoise_Subrange_np(ymatrix,smt_ymatrix,xwave,smt_xwave,startpnt,endpnt,suffix)
	wave ymatrix,smt_ymatrix,xwave,smt_xwave
	variable startpnt,endpnt
	string suffix
	
	variable nrows=dimsize(ymatrix,0)
	variable ncols=dimsize(ymatrix,1)
	variable nrows_smt=dimsize(smt_ymatrix,0)
	variable ncols_smt=dimsize(smt_ymatrix,1)
	
	if(ncols!=ncols_smt)
		print "Please make sure two matrices have the same number of columns"
		abort
	endif
	
	string ymatrixstr=nameofwave(ymatrix)
	string ymatrixsmtstr=nameofwave(smt_ymatrix)
	
	make/o/d/n=(nrows_smt,ncols_smt)/Free newymatrix
	make/o/d/n=(ncols_smt) $("AvgNoise"+suffix)
	make/o/d/n=(ncols_smt) $("NoiseSDV"+suffix)
	make/o/d/n=(nrows_smt)/Free tempwave

	wave newwave=$("AvgNoise"+suffix)
	wave newwave_sdv=$("NoiseSDV"+suffix)
	
	variable i,j
	
	for(i=0;i<ncols;i+=1)
		for(j=startpnt;j<endpnt;j+=1)
			findlevel/Q/P xwave,smt_xwave[j]
			newymatrix[j][i]=abs(ymatrix(V_levelx)[i]-smt_ymatrix[j][i])
		endfor
		tempwave=newymatrix[p][i]
		wavestats/q tempwave
		newwave[i]=V_avg
		newwave_sdv[i]=V_sdev
		
	endfor
	
	Note/k newwave "This wave stores the average difference between smoothed ( from matrix "+ymatrixstr+" ) and unsmoothed (from "+ymatrixsmtstr+" ) as noise level from row "+num2str(startpnt)+" to row "+num2str(endpnt)
	Note/k newwave_sdv "This wave stores the standard deviation of noise. Calculated along with the wave "+nameofwave(newwave)
	
End


//This function determines the threshold of noise and make a filter wave accordingly where points with noise>threshold =0; otherwise 1
Function FindNoiseThreshold()
	
	string noisewavestr = "AvgNoise_ar35"
	string masswavestr = "cmp_mass_rf_norm"
	variable Ionpercent = 5
	variable multiplier = 3
	string filterwavestr = "Filterwave_noise"
	
	prompt noisewavestr, "Select noise wave", popup WaveList("*",";","")
	prompt masswavestr, "Select mass contribution wave", popup WaveList("*",";","")
	prompt Ionpercent, "Set percent of ions with lowest noise"
	prompt multiplier, "Set a multiplier"
	prompt filterwavestr, "Name for the filterwave created"
	Doprompt "Select waves",noisewavestr, masswavestr, ionpercent, multiplier, filterwavestr
	
	if(V_flag)
		abort
	endif
	
	wave noisewave = $noisewavestr
	wave masswave = $masswavestr
	
	variable nrows=numpnts(noisewave)
	
	make/o/d/n=(nrows)/Free Noisewavesort 
	make/o/d/n=(nrows) $(filterwavestr)
	wave Filterwave_noise = $(filterwavestr)
	Filterwave_noise = 0
	
	Noisewavesort = noisewave
	
	sort noisewave,noisewavesort
	
	variable i_limit = ceil(nrows*Ionpercent/100)-1
	variable noise_limit = noisewavesort[i_limit]*multiplier
	
	//print i_limit
	
	variable i
	
	for(i=0;i<nrows;i+=1)
		if(noisewave[i] < noise_limit)
			filterwave_noise[i] = 1
		endif
	endfor
	
	make/o/d/n=(nrows)/Free masswave_fil
	masswave_fil = masswave
	
	masswave_fil = masswave[p] * filterwave_noise[p]
	wavestats/q masswave_fil
	
	print "The noise threshold is set to "+num2str(noise_limit)
	print "The remaining mass fraction after filtering out noisy ions is "+num2str(V_sum)
	
	Note/k filterwave_noise "The points in this wave = 1 if the corresponding values in wave "+noisewavestr+" < "+num2str(noise_limit)+" ; otherwise = 0"
	
End

Function FindNoiseThreshold_np(noisewave,masswave,Ionpercent,multiplier)
	wave noisewave, masswave
	variable Ionpercent, multiplier
	
	variable nrows=numpnts(masswave)
	
	make/o/d/n=(nrows)/Free Noisewavesort 
	make/o/d/n=(nrows) Filterwave_noise
	Filterwave_noise = 0
	
	Noisewavesort = noisewave
	
	sort noisewave,noisewavesort
	
	variable i_limit = ceil(nrows*Ionpercent/100)-1
	variable noise_limit = noisewavesort[i_limit]*multiplier
	
	//print i_limit
	
	variable i
	
	for(i=0;i<nrows;i+=1)
		if(noisewave[i] < noise_limit)
			filterwave_noise[i] = 1
		endif
	endfor
	
	make/o/d/n=(nrows)/Free masswave_fil
	masswave_fil = masswave
	
	masswave_fil = masswave[p] * filterwave_noise[p]
	wavestats/q masswave_fil
	
	print "The noise threshold is set to "+num2str(noise_limit)
	print "The remaining mass fraction after filtering out noisy ions is "+num2str(V_sum)
	
	Note/k filterwave_noise "The points in this wave = 1 if the corresponding values in wave "+nameofwave(noisewave)+" < "+num2str(noise_limit)+" ; otherwise = 0"
	
End

//This function create basic waves regarding the data set, including mass spectrum, bulk thermogram, m/z wave for each ion
Function Createwave_BasicP()
	
	string ymatrixstr = "cmp_matrix_rf"
	string nameswavestr = "Names_origin_rf"
	string suffix = "_rf"
	
	prompt ymatrixstr, "Select matrix of thermograms", popup WaveList("*",";","")
	prompt nameswavestr, "Select wave of compound names", popup WaveList("*",";","")
	prompt suffix, "Suffix to add to the new waves"
	Doprompt "Select waves",ymatrixstr, nameswavestr, suffix
	
	if(V_flag)
		abort
	endif
	
	wave ymatrix=$ymatrixstr
	wave/t nameswave=$nameswavestr
	
	variable nrows=dimsize(ymatrix,0)
	variable ncols=dimsize(ymatrix,1)
	
	make/o/d/n=(ncols) $("Molarmass_cmp"+suffix)=0
	make/o/d/n=(ncols) $("cmp_mass"+suffix)
	make/o/d/n=(ncols) $("cmp_mass_norm"+suffix)
	make/o/d/n=(ncols) $("cmp_mass_norm2Max"+suffix)
	make/o/d/n=(ncols) $("cmp_indx"+suffix)
	make/o/d/n=(nrows) $("Bulk_TG"+suffix)
	
	wave molarmasswave = $("Molarmass_cmp"+suffix)
	wave cmpmasswave = $("cmp_mass"+suffix)
	wave cmpmassnormwave = $("cmp_mass_norm"+suffix)
	wave cmpmassnorm2Mwave = $("cmp_mass_norm2Max"+suffix)
	wave cmpindxwave = $("cmp_indx"+suffix)
	wave bulkTGwave = $("Bulk_TG"+suffix)
	
	matrixop/o cmpmasswave = sumcols(ymatrix)^t
	matrixop/o bulkTGwave = sumrows(ymatrix)
	
	wavestats/q cmpmasswave
	cmpmassnormwave = cmpmasswave/V_sum
	cmpmassnorm2Mwave = cmpmasswave/V_max
	cmpindxwave = x
	
	string cmp_str
	variable i
	variable nC,nH,nO,nN,nS
	
	for(i=0;i<ncols;i+=1)
		cmp_str = nameswave[i]
		nC = ElementNumber(cmp_str,"C")
		nH = ElementNumber(cmp_str,"H")
		nO = ElementNumber(cmp_str,"O")
		nN = ElementNumber(cmp_str,"N")
		nS = ElementNumber(cmp_str,"S")
		molarmasswave[i]=12.000*nC+1.008*nH+15.995*nO+14.003*nN+31.972*nS
	endfor
	
	Note/k molarmasswave "This wave stores the molar mass of each ion according to name wave "+nameswavestr
	Note/k cmpmasswave "This wave stores the sum of each column of matrix "+ymatrixstr+" , corresponds to a mass spectrum"
	Note/k cmpmassnormwave "This wave is wave "+nameofwave(cmpmasswave)+" normalized to its sum"
	Note/k cmpmassnorm2Mwave "This wave is wave "+nameofwave(cmpmasswave)+" normalized to its max"
	Note/k cmpindxwave "This wave stores the index number of the ions in wave "+nameswavestr
	Note/k bulkTGwave "This wave stores the sum of each row of matrix "+ymatrixstr+" , corresponds to a thermogram of the bulk"
	
End

Function Createwave_BasicP_np(ymatrix, nameswave,suffix)
	
	wave ymatrix
	wave/t nameswave
	string suffix
	
	string nameswavestr = nameofwave(nameswave)
	string ymatrixstr = nameofwave(ymatrix)
	
	variable nrows=dimsize(ymatrix,0)
	variable ncols=dimsize(ymatrix,1)
	
	make/o/d/n=(ncols) $("Molarmass_cmp"+suffix)=0
	make/o/d/n=(ncols) $("cmp_mass"+suffix)
	make/o/d/n=(ncols) $("cmp_mass_norm"+suffix)
	make/o/d/n=(ncols) $("cmp_mass_norm2Max"+suffix)
	make/o/d/n=(ncols) $("cmp_indx"+suffix)
	make/o/d/n=(nrows) $("Bulk_TG"+suffix)
	
	wave molarmasswave = $("Molarmass_cmp"+suffix)
	wave cmpmasswave = $("cmp_mass"+suffix)
	wave cmpmassnormwave = $("cmp_mass_norm"+suffix)
	wave cmpmassnorm2Mwave = $("cmp_mass_norm2Max"+suffix)
	wave cmpindxwave = $("cmp_indx"+suffix)
	wave bulkTGwave = $("Bulk_TG"+suffix)
	
	matrixop/o cmpmasswave = sumcols(ymatrix)^t
	matrixop/o bulkTGwave = sumrows(ymatrix)
	
	wavestats/q cmpmasswave
	cmpmassnormwave = cmpmasswave/V_sum
	cmpmassnorm2Mwave = cmpmasswave/V_max
	cmpindxwave = x
	
	string cmp_str
	variable i
	variable nC,nH,nO,nN,nS
	
	for(i=0;i<ncols;i+=1)
		cmp_str = nameswave[i]
		nC = ElementNumber(cmp_str,"C")
		nH = ElementNumber(cmp_str,"H")
		nO = ElementNumber(cmp_str,"O")
		nN = ElementNumber(cmp_str,"N")
		nS = ElementNumber(cmp_str,"S")
		molarmasswave[i]=12.000*nC+1.008*nH+15.995*nO+14.003*nN+31.972*nS
	endfor
	
	Note/k molarmasswave "This wave stores the molar mass of each ion according to name wave "+nameswavestr
	Note/k cmpmasswave "This wave stores the sum of each column of matrix "+ymatrixstr+" , corresponds to a mass spectrum"
	Note/k cmpmassnormwave "This wave is wave "+nameofwave(cmpmasswave)+" normalized to its sum"
	Note/k cmpmassnorm2Mwave "This wave is wave "+nameofwave(cmpmasswave)+" normalized to its max"
	Note/k cmpindxwave "This wave stores the index number of the ions in wave "+nameswavestr
	Note/k bulkTGwave "This wave stores the sum of each row of matrix "+ymatrixstr+" , corresponds to a thermogram of the bulk"
	
End


//This function create basic waves regarding the data set, including mass spectrum, bulk thermogram, m/z wave, # of C, H, O, N, OSc, H/C, O/C for each ion
Function Createwave_BasicP_ext()
	
	string nameswavestr = "Names_fil"
	string suffix = "_fil"
	
	prompt nameswavestr, "Select wave of compound names", popup WaveList("*",";","")
	prompt suffix, "Suffix to add to the new waves"
	Doprompt "Select waves",nameswavestr, suffix
	
	if(V_flag)
		abort
	endif
	
	wave/t nameswave=$nameswavestr

	variable ncols=numpnts(nameswave)
	
	make/o/d/n=(ncols) $("Molarmass_cmp"+suffix)=0
	make/o/d/n=(ncols) $("C_num"+suffix)=0
	make/o/d/n=(ncols) $("H_num"+suffix)=0
	make/o/d/n=(ncols) $("O_num"+suffix)=0
	make/o/d/n=(ncols) $("N_num"+suffix)=0
	make/o/d/n=(ncols) $("S_num"+suffix)=0
	make/o/d/n=(ncols) $("HC"+suffix)=0
	make/o/d/n=(ncols) $("OC"+suffix)=0
	make/o/d/n=(ncols) $("NC"+suffix)=0
	make/o/d/n=(ncols) $("OSc"+suffix)=0
	
	wave molarmasswave = $("Molarmass_cmp"+suffix)
	wave nCwave = $("C_num"+suffix)
	wave nHwave = $("H_num"+suffix)
	wave nOwave = $("O_num"+suffix)
	wave nNwave = $("N_num"+suffix)
	wave nSwave = $("S_num"+suffix)
	wave HCrwave = $("HC"+suffix)
	wave OCrwave = $("OC"+suffix)
	wave NCrwave = $("NC"+suffix)
	wave OScwave = $("OSc"+suffix)
	
	string cmp_str
	variable i
	
	for(i=0;i<ncols;i+=1)
		cmp_str = nameswave[i]
		nCwave[i] = ElementNumber(cmp_str,"C")
		nHwave[i] = ElementNumber(cmp_str,"H")
		nOwave[i] = ElementNumber(cmp_str,"O")
		nNwave[i] = ElementNumber(cmp_str,"N")
		nSwave[i] = ElementNumber(cmp_str,"S")
		HCrwave[i] = nHwave[i]/nCwave[i]
		OCrwave[i] = nOwave[i]/nCwave[i]
		NCrwave[i] = nNwave[i]/nCwave[i]
		OScwave[i] = 2*OCrwave[i]-HCrwave[i]-5*NCrwave[i]
		molarmasswave[i]=12.000*nCwave[i]+1.008*nHwave[i]+15.995*nOwave[i]+14.003*nNwave[i]+31.972*nSwave[i]
	endfor
	
	Note/k nCwave "This waves stores the number of Carbon for each ion according to name wave "+nameswavestr
	Note/k nHwave "This waves stores the number of Hydrogen for each ion according to name wave "+nameswavestr
	Note/k nOwave "This waves stores the number of Oxygen for each ion according to name wave "+nameswavestr
	Note/k nNwave "This waves stores the number of Nitrogen for each ion according to name wave "+nameswavestr
	Note/k nSwave "This waves stores the number of Sulfur for each ion according to name wave "+nameswavestr
	Note/k HCrwave "This waves stores the H/C ratio for each ion according to name wave "+nameswavestr
	Note/k OCrwave "This waves stores the O/C ratio for each ion according to name wave "+nameswavestr
	Note/k NCrwave "This waves stores the N/C ratio for each ion according to name wave "+nameswavestr
	Note/k OScwave "This waves stores the oxidation state of Carbon for each ion according to name wave "+nameswavestr
	Note/k molarmasswave "This wave stores the molar mass of each ion according to name wave "+nameswavestr
	
End


Function Createwave_BasicP_ext_np(nameswave,suffix)
	
	wave/t nameswave
	string suffix
	
	string nameswavestr = nameofwave(nameswave)

	variable ncols=numpnts(nameswave)
	
	make/o/d/n=(ncols) $("Molarmass_cmp"+suffix)=0
	make/o/d/n=(ncols) $("C_num"+suffix)=0
	make/o/d/n=(ncols) $("H_num"+suffix)=0
	make/o/d/n=(ncols) $("O_num"+suffix)=0
	make/o/d/n=(ncols) $("N_num"+suffix)=0
	make/o/d/n=(ncols) $("S_num"+suffix)=0
	make/o/d/n=(ncols) $("HC"+suffix)=0
	make/o/d/n=(ncols) $("OC"+suffix)=0
	make/o/d/n=(ncols) $("NC"+suffix)=0
	make/o/d/n=(ncols) $("OSc"+suffix)=0

	
	wave molarmasswave = $("Molarmass_cmp"+suffix)
	wave nCwave = $("C_num"+suffix)
	wave nHwave = $("H_num"+suffix)
	wave nOwave = $("O_num"+suffix)
	wave nNwave = $("N_num"+suffix)
	wave nSwave = $("S_num"+suffix)
	wave HCrwave = $("HC"+suffix)
	wave OCrwave = $("OC"+suffix)
	wave NCrwave = $("NC"+suffix)
	wave OScwave = $("OSc"+suffix)
	
	string cmp_str
	variable i
	
	for(i=0;i<ncols;i+=1)
		cmp_str = nameswave[i]
		nCwave[i] = ElementNumber(cmp_str,"C")
		nHwave[i] = ElementNumber(cmp_str,"H")
		nOwave[i] = ElementNumber(cmp_str,"O")
		nNwave[i] = ElementNumber(cmp_str,"N")
		nSwave[i] = ElementNumber(cmp_str,"S")
		HCrwave[i] = nHwave[i]/nCwave[i]
		OCrwave[i] = nOwave[i]/nCwave[i]
		NCrwave[i] = nNwave[i]/nCwave[i]
		OScwave[i] = 2*OCrwave[i]-HCrwave[i]-5*NCrwave[i]
		molarmasswave[i]=12.000*nCwave[i]+1.008*nHwave[i]+15.995*nOwave[i]+14.003*nNwave[i]+31.972*nSwave[i]
	endfor
	
	Note/k nCwave "This waves stores the number of Carbon for each ion according to name wave "+nameswavestr
	Note/k nHwave "This waves stores the number of Hydrogen for each ion according to name wave "+nameswavestr
	Note/k nOwave "This waves stores the number of Oxygen for each ion according to name wave "+nameswavestr
	Note/k nNwave "This waves stores the number of Nitrogen for each ion according to name wave "+nameswavestr
	Note/k nSwave "This waves stores the number of Sulfur for each ion according to name wave "+nameswavestr
	Note/k HCrwave "This waves stores the H/C ratio for each ion according to name wave "+nameswavestr
	Note/k OCrwave "This waves stores the O/C ratio for each ion according to name wave "+nameswavestr
	Note/k NCrwave "This waves stores the N/C ratio for each ion according to name wave "+nameswavestr
	Note/k OScwave "This waves stores the oxidation state of Carbon for each ion according to name wave "+nameswavestr
	Note/k molarmasswave "This wave stores the molar mass of each ion according to name wave "+nameswavestr
		
End

// This function finds the number of a specific element in the compound
//c_str is the elemental composition of the compound
//e_str is the targeted element, e.g. C, H, O, I etc. 
//This function returns e_num which is the number of element e_str in the compound c_str
Function ElementNumber_str(c_str,e_str)
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


//This function determines the threshold for "significant" ions according the relationship between cumulative mass vs. individual mass
Function FindMassFraThreshold()

	string masswavestr = "cmp_mass_norm_rf"
	variable cumMassfra = 0.8
	
	prompt masswavestr, "Select wave of mass contribution", popup WaveList("*",";","")
	prompt cumMassfra, "Set threshold of cumulative mass fraction"
	Doprompt "Select waves",masswavestr,cumMassfra
	
	if(V_flag)
		abort
	endif
	
	wave masswave = $masswavestr
	
	variable pntN=numpnts(masswave)
	
	make/o/d/n=(pntN) $(nameofwave(masswave)+"_sort")
	make/o/d/n=(pntN) $(nameofwave(masswave)+"_INT")
	wave masswavesort=$(nameofwave(masswave)+"_sort")
	wave masswaveint=$(nameofwave(masswave)+"_INT")
	masswavesort=masswave
	
	string notestr = "Function FindMassFraThreshold is used. Masswave <- "+nameofwave(masswave)+" and cumMassfra = "+num2str(cumMassfra)
	
	sort/r masswave,masswavesort
	Integrate masswavesort/D=masswaveint
	
	findlevel/q/p masswaveint,cumMassfra
	variable threshold=masswavesort[V_levelx]
	
	extract/indx/o masswave, HM_indx,masswave>=threshold
	
	Note/k HM_indx "This wave stores the index of the ions that have mass fraction larger than "+num2str(threshold)
	Note HM_indx notestr
	Note/k masswavesort "This wave stores the sorted values of "+nameofwave(masswave)+" from large to small"
	Note masswavesort notestr
	Note/k masswaveint "This wave stores the cumulative mass fraction according to wave "+nameofwave(masswavesort)
	Note masswaveint notestr
	
	print "The threshold for high-mass ions is: "+num2str(threshold)
	print "The number of high-mass ions is: "+num2str(numpnts(hm_indx))
	return threshold
	
End

Function FindMassFraThreshold_np(masswave,cumMassfra)
	wave masswave
	variable cumMassfra
	
	variable pntN=numpnts(masswave)
	
	make/o/d/n=(pntN) $(nameofwave(masswave)+"_sort")
	make/o/d/n=(pntN) $(nameofwave(masswave)+"_INT")
	wave masswavesort=$(nameofwave(masswave)+"_sort")
	wave masswaveint=$(nameofwave(masswave)+"_INT")
	masswavesort=masswave
	
	string notestr = "Function FindMassFraThreshold_np is used. Masswave <- "+nameofwave(masswave)+" and cumMassfra = "+num2str(cumMassfra)
	
	sort/r masswave,masswavesort
	Integrate masswavesort/D=masswaveint
	
	findlevel/q/p masswaveint,cumMassfra
	variable threshold=masswavesort[V_levelx]
	
	extract/indx/o masswave, HM_indx,masswave>=threshold
	
	Note/k HM_indx "This wave stores the index of the ions that have mass fraction larger than "+num2str(threshold)
	Note HM_indx notestr
	Note/k masswavesort "This wave stores the sorted wave of "+nameofwave(masswave)+" from large to small"
	Note masswavesort notestr
	Note/k masswaveint "This wave stores the cumulative mass fraction according to wave "+nameofwave(masswavesort)
	Note masswaveint notestr
	
	return threshold
	
	
End


//This function filters out columns in matrix when corresponding column index# in filterwve is 0
//Length of the filter should be the same as number of columns of the matrix
Function  FilterMatrix_Col()
	
	string ywavestr = "cmp_matrix_rf_smtnorm_10a"
	string filterstr = "Filterwave_noise"
	string newwavestr = "cmp_matrix_sn_10a_fil"
	
	prompt ywavestr, "Select matrix to filter", popup WaveList("*",";","")
	prompt filterstr, "Select filter wave", popup WaveList("*",";","")
	prompt newwavestr, "Name of the filtered matrix"
	Doprompt "Select waves",ywavestr,filterstr,newwavestr
	
	if(V_flag)
		abort
	endif
	
	wave ywave = $ywavestr
	wave filter = $filterstr
	
	variable nrows=dimsize(ywave,0)
	variable ncols=dimsize(ywave,1)
	
	variable i
	variable count
	count=0
	
	for(i=0;i<ncols;i+=1)
		if(filter[i]==1)
			count+=1
		endif
	endfor
	
	make/o/d/n=(nrows,count) $newwavestr
	wave newwave=$newwavestr
	
	count=0
	for(i=0;i<ncols;i+=1)
		if(filter[i]==1)
			newwave[][count]=ywave[p][i]
			count+=1
		endif
	endfor
	
	string notestr
	notestr="This matrix stores the columns in the wave "+nameofwave(ywave)+" when corresponding values in the wave "+nameofwave(filter)+" = 1"
	Note/k newwave notestr

End


Function  FilterMatrix_Col_np(ywave,filter,newwavestr)
	
	wave ywave
	wave filter
	string newwavestr
	
	variable nrows=dimsize(ywave,0)
	variable ncols=dimsize(ywave,1)
	
	variable i
	variable count
	count=0
	
	for(i=0;i<ncols;i+=1)
		if(filter[i]==1)
			count+=1
		endif
	endfor
	
	make/o/d/n=(nrows,count) $newwavestr
	wave newwave=$newwavestr
	
	count=0
	for(i=0;i<ncols;i+=1)
		if(filter[i]==1)
			newwave[][count]=ywave[p][i]
			count+=1
		endif
	endfor
	
	string notestr
	notestr="This matrix stores the columns in the wave "+nameofwave(ywave)+" when corresponding values in the wave "+nameofwave(filter)+" = 1"
	Note/k newwave notestr

End

//This function filters out points in wave when corresponding index in filterwve is 0
//The wave is a numeric wave
//See FilterTextWave for filtering text wave based on filterwave
Function  Filter1DWave()
	
	string ywavestr = "AvgNoise_ar35"
	string filterstr = "Filterwave_noise"
	string newwavestr ="AvgNoise_ar35_fil"
	
	prompt ywavestr, "Select numeric wave to filter", popup WaveList("*",";","")
	prompt filterstr, "Select filter wave", popup WaveList("*",";","")
	prompt newwavestr, "Name of the filtered wave"
	Doprompt "Select waves",ywavestr,filterstr,newwavestr
	
	if(V_flag)
		abort
	endif
	
	wave ywave = $ywavestr
	wave filter = $filterstr
	
	variable nrows=numpnts(filter)
	
	variable i
	variable count
	count=0
	
	for(i=0;i<nrows;i+=1)
		if(filter[i]==1)
			count+=1
		endif
	endfor
	
	make/o/d/n=(count) $newwavestr
	wave newwave=$newwavestr
	
	count=0
	for(i=0;i<nrows;i+=1)
		if(filter[i]==1)
			newwave[count]=ywave[i]
			count+=1
		endif
	endfor
	
	string notestr
	notestr="This wave stores values in the wave "+ywavestr+" when corresponding values in the wave "+filterstr+" = 1"
	Note/k newwave notestr

End

Function  Filter1DWave_np(ywave,filter,newwavestr)
	
	wave ywave
	wave filter
	string newwavestr
	
	variable nrows=numpnts(filter)
	
	variable i
	variable count
	count=0
	
	for(i=0;i<nrows;i+=1)
		if(filter[i]==1)
			count+=1
		endif
	endfor
	
	make/o/d/n=(count) $newwavestr
	wave newwave=$newwavestr
	
	count=0
	for(i=0;i<nrows;i+=1)
		if(filter[i]==1)
			newwave[count]=ywave[i]
			count+=1
		endif
	endfor
	
	string notestr
	notestr="This wave stores values in the wave "+nameofwave(ywave)+" when corresponding values in the wave "+nameofwave(filter)+" = 1"
	Note/k newwave notestr

End


//This function filters out points in wave when corresponding index in filterwve is 0
//The wave is a text wave
Function  FilterTextWave()
	
	string ywavestr = "Names_origin_rf"
	string filterstr = "Filterwave_noise"
	string newwavestr ="Names_origin_fil"
	
	prompt ywavestr, "Select text wave to filter", popup WaveList("*",";","")
	prompt filterstr, "Select filter wave", popup WaveList("*",";","")
	prompt newwavestr, "Name of the filtered wave"
	Doprompt "Select waves",ywavestr,filterstr,newwavestr
	
	if(V_flag)
		abort
	endif
	
	wave/t ywave = $ywavestr
	wave filter = $filterstr
	
	variable nrows=numpnts(filter)
	
	variable i
	variable count
	count=0
	
	for(i=0;i<nrows;i+=1)
		if(filter[i]==1)
			count+=1
		endif
	endfor
	
	make/o/t/n=(count) $newwavestr
	wave/t newwave=$newwavestr
	
	count=0
	for(i=0;i<nrows;i+=1)
		if(filter[i]==1)
			newwave[count]=ywave[i]
			count+=1
		endif
	endfor
	
	string notestr
	notestr="This wave stores values in the wave "+ywavestr+" when corresponding values in the wave "+filterstr+" = 1"
	Note/k newwave notestr

End


Function  FilterTextWave_np(ywave,filter,newwavestr)
	
	wave/t ywave
	wave filter
	string newwavestr
	
	variable nrows=numpnts(filter)
	
	variable i
	variable count
	count=0
	
	for(i=0;i<nrows;i+=1)
		if(filter[i]==1)
			count+=1
		endif
	endfor
	
	make/o/t/n=(count) $newwavestr
	wave/t newwave=$newwavestr
	
	count=0
	for(i=0;i<nrows;i+=1)
		if(filter[i]==1)
			newwave[count]=ywave[i]
			count+=1
		endif
	endfor
	
	string notestr
	notestr="This wave stores values in the wave "+nameofwave(ywave)+" when corresponding values in the wave "+nameofwave(filter)+" = 1"
	Note/k newwave notestr

End

//This function calculates the Euclidean distance all the column pairs in ymatrix
Function CalculateED()

	string ymatrixstr = "cmp_matrix_sn_10a_fil"
	string newEDwavestr = "ED_fil"
	
	prompt ymatrixstr, "Select matrix", popup WaveList("*",";","")
	prompt newEDwavestr, "Name of the new matrix with ED"
	Doprompt "Select waves",ymatrixstr,newEDwavestr
	
	if(V_flag)
		abort
	endif

	wave ymatrix = $ymatrixstr
	
	variable nrows=dimsize(ymatrix,0)
	variable ncols=dimsize(ymatrix,1)
	variable i,j
	
	make/o/d/n=(nrows)/Free this,that,squared
	make/o/d/n=(ncols,ncols) $(newEDwavestr)
	
	wave EuclideanDwave=$(newEDwavestr)
	
	for(i=0;i<ncols;i+=1)
		this=ymatrix[p][i]
		for(j=0;j<ncols;j+=1)
			that=ymatrix[p][j]
			squared=this-that
			squared=squared^2
			wavestats/q squared
			EuclideanDwave[i][j]=V_sum^0.5
		endfor
	endfor
	
	Note/k EuclideanDwave "This matrix stores the Euclidean distance between all the pairs of columns in matrix "+nameofwave(ymatrix)
End


Function CalculateED_np(ymatrix,newEDwavestr)
	wave ymatrix
	string newEDwavestr
	
	variable nrows=dimsize(ymatrix,0)
	variable ncols=dimsize(ymatrix,1)
	variable i,j
	
	make/o/d/n=(nrows)/Free this,that,squared
	make/o/d/n=(ncols,ncols) $(newEDwavestr)
	
	wave EuclideanDwave=$(newEDwavestr)
	
	for(i=0;i<ncols;i+=1)
		this=ymatrix[p][i]
		for(j=0;j<ncols;j+=1)
			that=ymatrix[p][j]
			squared=this-that
			squared=squared^2
			wavestats/q squared
			EuclideanDwave[i][j]=V_sum^0.5
		endfor
	endfor
	
	Note/k EuclideanDwave "This matrix stores the Euclidean distance between all the pairs of columns in matrix "+nameofwave(ymatrix)
End


//This function executesthe NSSClustering without the consideration of unclustered high-mass ions as additional one-member clusters
//This function use a specified epsilon
Function NSSC()

	string matstr = "cmp_matrix_sn_10a_fil"
	string Dmatstr = "ED_fil"
	string sortingwavestr = "AvgNoise_ar35_fil"
	string weightwavestr = "cmp_mass_norm_fil"
	variable nclusters = 20
	variable eps = 2.5
	
	prompt matstr, "Select matrix of thermograms", popup WaveList("*",";","")
	prompt Dmatstr, "Select matrix of ED", popup WaveList("*",";","")
	prompt sortingwavestr, "Select noise wave", popup WaveList("*",";","")
	prompt weightwavestr, "Select mass contribution wave", popup WaveList("*",";","")
	prompt nclusters, "Max # of clusters expected"
	prompt eps, "Enter distance crriterion eps"
	Doprompt "Select waves",matstr, Dmatstr, sortingwavestr, weightwavestr, nclusters, eps
	
	if(V_flag)
		abort
	endif
	
	wave mat = $matstr
	wave Dmat = $Dmatstr
	wave sortingwave = $sortingwavestr
	wave weightwave = $weightwavestr
	
	
	variable nrows = dimsize(mat,0) // length
	variable ncols = dimsize(mat,1) // each spectra
	
	variable minneighbors=2 // the minimum number of neighbors for something to have to be considered a cluster

	
	variable i,j, k
	
	make/o/d/n=2/Free ungroupedN
	
	make/o/d/n=(ncols)/Free Pstate = 0 // CLUSTER state
	make/o/d/n=(ncols)/Free Vstate = 0 // SCAN state
	make/o/d/n=(ncols) Clusternum=-1
	make/o/d/n=(nrows,nclusters) Clusterseed=0
	make/o/d/n=(nclusters) ClusterseedIndx=0
	
	make/o/d/n=(nrows)/Free this,that,squared
	duplicate/o mat mat_temp
	
	variable cc = 0
	variable dotvar, dotref
	variable seed_idx
	make/o/d/n=(ncols)/Free DotArray = nan
	
	k = 0
	cc = 0
	do
		cc+=1
		if(cc > 2*ncols)
			//print "too many loops. Probably means not enough neighbors."
			break
		endif
		extract/INDX/o vstate, ptemp, vstate==0 // extract indices of all spectra that have not been checked yet
		extract/o sortingwave, sorttemp,vstate==0
		if(numpnts(ptemp)==0)
			break
		endif
		wavestats/q sorttemp
		seed_idx=V_minloc  //Used when sorted by noise, low->high
		//seed_idx=V_maxloc //Used when sorted by mass contribution, high->low
		//seed_idx = round(abs(enoise(numpnts(ptemp)-1))) // determine a random seed from this list
		seed_idx = ptemp[seed_idx] // get the index number, specifically, of the random seed
		vstate[seed_idx]=1
		dotref = dmat[seed_idx][seed_idx] // the reference dot product
		DotArray = dmat[p][seed_idx] // the dot products for the single spectrum to compare against
		//pstate[seed_idx] = 1 // set pstate to 1 for this seed #
		// check number of neighbors by extracting all that are within range eps
		extract/INDX/o DotArray, neighbors, abs(dotarray - dotref) < eps && pstate==0 
		
		if(numpnts(neighbors) < minNeighbors)
			continue // if there is no match, return to start
		endif
		
		
		cc +=1
		
		for(i=0;i<numpnts(neighbors);i+=1)
			pstate[neighbors[i]] = 1
			vstate[neighbors[i]] = 1
			Clusternum[neighbors[i]]=k+1
		endfor

		Clusterseed[][k]=mat[p][seed_idx]
		ClusterseedIndx[k]=seed_idx
			
		k+=1
	while(k<nclusters)

	extract/INDX/o pstate, ptemp, pstate==0
	//print numpnts(ptemp)
	ungroupedN[0]=numpnts(ptemp)
	mat=mat_temp	
	
	killwaves neighbors,mat_temp,ptemp,sorttemp
	
	duplicate/o clusternum Clusternum_org	
	
	variable reclusternum = 1
	variable unclusternum = 0
	
	j=0
	
	do
		ClusteronClusterNum_Wt_np(Clusternum,mat,weightwave,"_wt")
		wave Clusterave_wt
	
		extract/Indx/o clusternum,unclustered_indx,clusternum==-1
		unclusternum = numpnts(unclustered_indx)
		
		
		if(unclusternum != numpnts(clusternum))
			make/o/d/n=(unclusternum) ClosestClstED = 0
			make/o/d/n=(unclusternum) ClosestClstindx = 0
		
			for(j=0; j<unclusternum;j+=1)
				FindED_ClosestCluster(unclustered_indx[j],mat,Clusterave_wt)
				wave EDindx
				ClosestClstED[j] = EDindx[1]
				ClosestClstindx[j] = EDindx[0]
			endfor
		
			extract/indx/o closestclstED,ReCluster_indx,ClosestClstED < eps
			extract/o closestclstindx,ClosestClstindx_eps,ClosestClstED<eps
	
			reclusternum=numpnts(Recluster_indx)
	
	
			for(i=0;i<reclusternum;i+=1)
				clusternum[unclustered_indx[Recluster_indx[i]]]=ClosestClstindx_eps[i]+1
			endfor
		else
			reclusternum = 0
		endif
		
	while(reclusternum!=0) 
	
	wave ClusterSdv_wt,clusterupper_wt,clusterlower_wt
	killwaves recluster_indx,closestclstindx_eps, ClusterSdv_wt,clusterupper_wt,clusterlower_wt
	killwaves unclustered_indx,EDindx,ClosestclstED, ClosestclstIndx
	
	string Notestr="Function used: NSSC; mat <- "+nameofwave(mat)+" ; Dmat <- "+nameofwave(Dmat)+" ; sortingwave <- "+nameofwave(sortingwave)+" ; Weightwave<- "+nameofwave(Weightwave)+" ; nclusters = "+num2istr(nclusters)+" ; eps = "+num2str(eps)
	
	Note/k clusternum "This wave stores Cluster# for each ion, with -1 indicating unclustered"
	Note clusternum Notestr
	Note/k clusterseedindx "This wave stores the index of seed for each cluster"
	Note clusterseedindx Notestr
	Note/k clusternum_org "This wave stores Cluster# for each ion before executing the cycles to recluster unclustered ions"
	Note clusternum_org Notestr
	Note/k Clusterseed "Each column of this matrix stores the thermogram of the seed used for clustering. Column 0 corresponds to Cluster1"
	Note clusterseed Notestr
	
end

Function NSSC_np(mat,Dmat,sortingwave,weightwave, nclusters,eps)
	wave mat //cmp_matrix
	wave Dmat  //ED_matrix
	wave sortingwave
	wave weightwave
	variable nclusters  //the maximum number of clusters you would expect. 
	variable eps // 0.1 // the tolerance on the dot product. 
	
	
	variable nrows = dimsize(mat,0) // length
	variable ncols = dimsize(mat,1) // each spectra
	
	variable minneighbors=2 // the minimum number of neighbors for something to have to be considered a cluster

	
	variable i,j, k
	
	make/o/d/n=2/Free ungroupedN
	
	make/o/d/n=(ncols)/Free Pstate = 0 // CLUSTER state
	make/o/d/n=(ncols)/Free Vstate = 0 // SCAN state
	make/o/d/n=(ncols) Clusternum=-1
	make/o/d/n=(nrows,nclusters) Clusterseed=0
	make/o/d/n=(nclusters) ClusterseedIndx=0
	
	make/o/d/n=(nrows)/Free this,that,squared
	duplicate/o mat mat_temp
	
	variable cc = 0
	variable dotvar, dotref
	variable seed_idx
	make/o/d/n=(ncols)/Free DotArray = nan
	
	k = 0
	cc = 0
	do
		cc+=1
		if(cc > 2*ncols)
			//print "too many loops. Probably means not enough neighbors."
			break
		endif
		extract/INDX/o vstate, ptemp, vstate==0 // extract indices of all spectra that have not been checked yet
		extract/o sortingwave, sorttemp,vstate==0
		if(numpnts(ptemp)==0)
			break
		endif
		wavestats/q sorttemp
		seed_idx=V_minloc  //Used when sorted by noise, low->high
		//seed_idx=V_maxloc //Used when sorted by mass contribution, high->low
		//seed_idx = round(abs(enoise(numpnts(ptemp)-1))) // determine a random seed from this list
		seed_idx = ptemp[seed_idx] // get the index number, specifically, of the random seed
		vstate[seed_idx]=1
		dotref = dmat[seed_idx][seed_idx] // the reference dot product
		DotArray = dmat[p][seed_idx] // the dot products for the single spectrum to compare against
		//pstate[seed_idx] = 1 // set pstate to 1 for this seed #
		// check number of neighbors by extracting all that are within range eps
		extract/INDX/o DotArray, neighbors, abs(dotarray - dotref) < eps && pstate==0 
		
		if(numpnts(neighbors) < minNeighbors)
			continue // if there is no match, return to start
		endif
		
		
		cc +=1
		
		for(i=0;i<numpnts(neighbors);i+=1)
			pstate[neighbors[i]] = 1
			vstate[neighbors[i]] = 1
			Clusternum[neighbors[i]]=k+1
		endfor

		Clusterseed[][k]=mat[p][seed_idx]
		ClusterseedIndx[k]=seed_idx
			
		k+=1
	while(k<nclusters)

	extract/INDX/o pstate, ptemp, pstate==0
	//print numpnts(ptemp)
	ungroupedN[0]=numpnts(ptemp)
	mat=mat_temp	
	
	killwaves neighbors,mat_temp,ptemp,sorttemp
	
	duplicate/o clusternum Clusternum_org	
	
	variable reclusternum = 1
	variable unclusternum = 0
	
	j=0
	
	do
		ClusteronClusterNum_Wt_np(Clusternum,mat,weightwave,"_wt")
		wave Clusterave_wt
	
		extract/Indx/o clusternum,unclustered_indx,clusternum==-1
		unclusternum = numpnts(unclustered_indx)
		
		if(unclusternum != numpnts(clusternum))
			make/o/d/n=(unclusternum) ClosestClstED = 0
			make/o/d/n=(unclusternum) ClosestClstindx = 0
		
			for(j=0; j<unclusternum;j+=1)
				FindED_ClosestCluster(unclustered_indx[j],mat,Clusterave_wt)
				wave EDindx
				ClosestClstED[j] = EDindx[1]
				ClosestClstindx[j] = EDindx[0]
			endfor
		
			extract/indx/o closestclstED,ReCluster_indx,ClosestClstED < eps
			extract/o closestclstindx,ClosestClstindx_eps,ClosestClstED<eps
	
			reclusternum=numpnts(Recluster_indx)
	
	
			for(i=0;i<reclusternum;i+=1)
				clusternum[unclustered_indx[Recluster_indx[i]]]=ClosestClstindx_eps[i]+1
			endfor
		else
			reclusternum = 0
		endif
		
		//j+=1
	while(reclusternum!=0) 
	
	//print j
	killwaves recluster_indx,closestclstindx_eps
	killwaves unclustered_indx,EDindx,ClosestclstED, ClosestclstIndx
	
	string Notestr="Function used: NSSC_np; mat <- "+nameofwave(mat)+" ; Dmat <- "+nameofwave(Dmat)+" ; sortingwave <- "+nameofwave(sortingwave)+" ; Weightwave<- "+nameofwave(Weightwave)+" ; nclusters = "+num2istr(nclusters)+" ; eps = "+num2str(eps)
	
	Note/k clusternum "This wave stores Cluster# for each ion, with -1 indicating unclustered"
	Note clusternum Notestr
	Note/k clusterseedindx "This wave stores the index of seed for each cluster"
	Note clusterseedindx Notestr
	Note/k clusternum_org "This wave stores Cluster# for each ion before executing the cycles to recluster unclustered ions"
	Note clusternum_org Notestr
	Note/k Clusterseed "Each column of this matrix stores the thermogram of the seed used for clustering. Column 0 corresponds to Cluster1"
	Note clusterseed Notestr
	
end

//This function extends function NSSC to include one-member clusters that are significant high-mass-contribution unclustered ions
Function NSSC_w1M()

	string matstr = "cmp_matrix_sn_10a_fil"
	string Dmatstr = "ED_fil"
	string sortingwavestr = "AvgNoise_ar35_fil"
	string weightwavestr = "cmp_mass_norm_fil"
	variable nclusters = 50
	variable eps = 2.5
	variable Mthreshold = 0.05
	
	prompt matstr, "Select matrix of thermograms", popup WaveList("*",";","")
	prompt Dmatstr, "Select matrix of ED", popup WaveList("*",";","")
	prompt sortingwavestr, "Select noise wave", popup WaveList("*",";","")
	prompt weightwavestr, "Select mass contribution wave", popup WaveList("*",";","")
	prompt nclusters, "Max # of clusters expected"
	prompt eps, "Enter distance crriterion eps"
	prompt Mthreshold, "Threshold for significant ions"
	Doprompt "Select waves",matstr, Dmatstr, sortingwavestr, weightwavestr, nclusters, eps, Mthreshold
	
	if(V_flag)
		abort
	endif
	
	wave mat = $matstr
	wave Dmat = $Dmatstr
	wave sortingwave = $sortingwavestr
	wave weightwave = $weightwavestr
	
	
	variable nrows = dimsize(mat,0) // length
	variable ncols = dimsize(mat,1) // each spectra
	
	variable minneighbors=2 // the minimum number of neighbors for something to have to be considered a cluster

	
	variable i,j, k
	
	make/o/d/n=2/Free ungroupedN
	
	make/o/d/n=(ncols)/Free Pstate = 0 // CLUSTER state
	make/o/d/n=(ncols)/Free Vstate = 0 // SCAN state
	make/o/d/n=(ncols) Clusternum=-1
	make/o/d/n=(nrows,nclusters) Clusterseed=0
	make/o/d/n=(nclusters) ClusterseedIndx=0
	
	make/o/d/n=(nrows)/Free this,that,squared
	duplicate/o mat mat_temp
	
	variable cc = 0
	variable dotvar, dotref
	variable seed_idx
	make/o/d/n=(ncols)/Free DotArray = nan
	
	k = 0
	cc = 0
	do
		cc+=1
		if(cc > 2*ncols)
			//print "too many loops. Probably means not enough neighbors."
			break
		endif
		extract/INDX/o vstate, ptemp, vstate==0 // extract indices of all spectra that have not been checked yet
		extract/o sortingwave, sorttemp,vstate==0
		if(numpnts(ptemp)==0)
			break
		endif
		wavestats/q sorttemp
		seed_idx=V_minloc  //Used when sorted by noise, low->high
		//seed_idx=V_maxloc //Used when sorted by mass contribution, high->low
		//seed_idx = round(abs(enoise(numpnts(ptemp)-1))) // determine a random seed from this list
		seed_idx = ptemp[seed_idx] // get the index number, specifically, of the random seed
		vstate[seed_idx]=1
		dotref = dmat[seed_idx][seed_idx] // the reference dot product
		DotArray = dmat[p][seed_idx] // the dot products for the single spectrum to compare against
		//pstate[seed_idx] = 1 // set pstate to 1 for this seed #
		// check number of neighbors by extracting all that are within range eps
		extract/INDX/o DotArray, neighbors, abs(dotarray - dotref) < eps && pstate==0 
		
		if(numpnts(neighbors) < minNeighbors)
			continue // if there is no match, return to start
		endif
		
		
		cc +=1
		
		for(i=0;i<numpnts(neighbors);i+=1)
			pstate[neighbors[i]] = 1
			vstate[neighbors[i]] = 1
			Clusternum[neighbors[i]]=k+1
		endfor

		Clusterseed[][k]=mat[p][seed_idx]
		ClusterseedIndx[k]=seed_idx
			
		k+=1
	while(k<nclusters)

	extract/INDX/o pstate, ptemp, pstate==0
	//print numpnts(ptemp)
	ungroupedN[0]=numpnts(ptemp)
	mat=mat_temp	
	
	killwaves neighbors,mat_temp,ptemp,sorttemp
	
	duplicate/o clusternum Clusternum_org	
	
	variable reclusternum = 1
	variable unclusternum = 0
	
	j=0
	
	do
		ClusteronClusterNum_Wt_np(Clusternum,mat,weightwave,"_wt")
		wave Clusterave_wt
	
		extract/Indx/o clusternum,unclustered_indx,clusternum==-1
		unclusternum = numpnts(unclustered_indx)
		
		make/o/d/n=(unclusternum) ClosestClstED = 0
		make/o/d/n=(unclusternum) ClosestClstindx = 0
		
		for(j=0; j<unclusternum;j+=1)
			FindED_ClosestCluster(unclustered_indx[j],mat,Clusterave_wt)
			wave EDindx
			ClosestClstED[j] = EDindx[1]
			ClosestClstindx[j] = EDindx[0]
		endfor
		
		extract/indx/o closestclstED,ReCluster_indx,ClosestClstED < eps
		extract/o closestclstindx,ClosestClstindx_eps,ClosestClstED<eps
	
		reclusternum=numpnts(Recluster_indx)
	
	
		for(i=0;i<reclusternum;i+=1)
			clusternum[unclustered_indx[Recluster_indx[i]]]=ClosestClstindx_eps[i]+1
		endfor
		
		j+=1
	while(reclusternum!=0) 
	
	extract/indx/o clusternum, unclustered_hm_indx,clusternum==-1 && weightwave >= Mthreshold 
	
	variable nunclstHM = numpnts(unclustered_hm_indx)
	
	wavestats/q clusternum
	variable ncurrentclst = V_max
	variable newclstnum = ncurrentclst+1
	
	for(j=0 ; j<nunclstHM ;j+=1)
		clusternum[unclustered_hm_indx[j]] = newclstnum
		if(newclstnum-1 >= nclusters)
			Insertpoints newclstnum-1,1,clusterseedindx
			Insertpoints/M=1 newclstnum-1,1,clusterseed
		endif
		clusterseedindx[newclstnum-1] = unclustered_hm_indx[j]
		clusterseed[][newclstnum-1] = mat[p][unclustered_hm_indx[j]]
		newclstnum+=1
	endfor
	
	ClusteronClusterNum_Wt_np(Clusternum,mat,weightwave,"_wt")
	
	wave ClusterSdv_wt,clusterupper_wt,clusterlower_wt
	killwaves recluster_indx,closestclstindx_eps, ClusterSdv_wt,clusterupper_wt,clusterlower_wt
	killwaves unclustered_indx,unclustered_hm_indx,EDindx,ClosestclstED, ClosestclstIndx
	
	string Notestr="Function used: NSSC_w1M; mat <- "+nameofwave(mat)+" ; Dmat <- "+nameofwave(Dmat)+" ; sortingwave <- "+nameofwave(sortingwave)+" ; Weightwave<- "+nameofwave(Weightwave)+" ; nclusters = "+num2istr(nclusters)+" ; eps = "+num2str(eps)
	
	Note/k clusternum "This wave stores Cluster# for each ion, with -1 indicating unclustered"
	Note clusternum Notestr
	Note/k clusterseedindx "This wave stores the index of seed for each cluster"
	Note clusterseedindx Notestr
	Note/k clusternum_org "This wave stores Cluster# for each ion before executing the cycles to recluster unclustered ions"
	Note clusternum_org Notestr
	Note/k Clusterseed "Each column of this matrix stores the thermogram of the seed used for clustering. Column 0 corresponds to Cluster1"
	Note clusterseed Notestr
	
end

Function NSSC_w1M_np(mat,Dmat,sortingwave,weightwave,nclusters,eps,Mthreshold)

	wave mat, Dmat, sortingwave, weightwave
	variable nclusters, eps,Mthreshold
	
	variable nrows = dimsize(mat,0) // length
	variable ncols = dimsize(mat,1) // each spectra
	
	variable minneighbors=2 // the minimum number of neighbors for something to have to be considered a cluster

	
	variable i,j, k
	
	make/o/d/n=2/Free ungroupedN
	
	make/o/d/n=(ncols)/Free Pstate = 0 // CLUSTER state
	make/o/d/n=(ncols)/Free Vstate = 0 // SCAN state
	make/o/d/n=(ncols) Clusternum=-1
	make/o/d/n=(nrows,nclusters) Clusterseed=0
	make/o/d/n=(nclusters) ClusterseedIndx=0
	
	make/o/d/n=(nrows)/Free this,that,squared
	duplicate/o mat mat_temp
	
	variable cc = 0
	variable dotvar, dotref
	variable seed_idx
	make/o/d/n=(ncols)/Free DotArray = nan
	
	k = 0
	cc = 0
	do
		cc+=1
		if(cc > 2*ncols)
			//print "too many loops. Probably means not enough neighbors."
			break
		endif
		extract/INDX/o vstate, ptemp, vstate==0 // extract indices of all spectra that have not been checked yet
		extract/o sortingwave, sorttemp,vstate==0
		if(numpnts(ptemp)==0)
			break
		endif
		wavestats/q sorttemp
		seed_idx=V_minloc  //Used when sorted by noise, low->high
		//seed_idx=V_maxloc //Used when sorted by mass contribution, high->low
		//seed_idx = round(abs(enoise(numpnts(ptemp)-1))) // determine a random seed from this list
		seed_idx = ptemp[seed_idx] // get the index number, specifically, of the random seed
		vstate[seed_idx]=1
		dotref = dmat[seed_idx][seed_idx] // the reference dot product
		DotArray = dmat[p][seed_idx] // the dot products for the single spectrum to compare against
		//pstate[seed_idx] = 1 // set pstate to 1 for this seed #
		// check number of neighbors by extracting all that are within range eps
		extract/INDX/o DotArray, neighbors, abs(dotarray - dotref) < eps && pstate==0 
		
		if(numpnts(neighbors) < minNeighbors)
			continue // if there is no match, return to start
		endif
		
		
		cc +=1
		
		for(i=0;i<numpnts(neighbors);i+=1)
			pstate[neighbors[i]] = 1
			vstate[neighbors[i]] = 1
			Clusternum[neighbors[i]]=k+1
		endfor

		Clusterseed[][k]=mat[p][seed_idx]
		ClusterseedIndx[k]=seed_idx
			
		k+=1
	while(k<nclusters)

	extract/INDX/o pstate, ptemp, pstate==0
	//print numpnts(ptemp)
	ungroupedN[0]=numpnts(ptemp)
	mat=mat_temp	
	
	killwaves neighbors,mat_temp,ptemp,sorttemp
	
	duplicate/o clusternum Clusternum_org	
	
	variable reclusternum = 1
	variable unclusternum = 0
	
	j=0
	
	do
		ClusteronClusterNum_Wt_np(Clusternum,mat,weightwave,"_wt")
		wave Clusterave_wt
	
		extract/Indx/o clusternum,unclustered_indx,clusternum==-1
		unclusternum = numpnts(unclustered_indx)
		
		make/o/d/n=(unclusternum) ClosestClstED = 0
		make/o/d/n=(unclusternum) ClosestClstindx = 0
		
		for(j=0; j<unclusternum;j+=1)
			FindED_ClosestCluster(unclustered_indx[j],mat,Clusterave_wt)
			wave EDindx
			ClosestClstED[j] = EDindx[1]
			ClosestClstindx[j] = EDindx[0]
		endfor
		
		extract/indx/o closestclstED,ReCluster_indx,ClosestClstED < eps
		extract/o closestclstindx,ClosestClstindx_eps,ClosestClstED<eps
	
		reclusternum=numpnts(Recluster_indx)
	
	
		for(i=0;i<reclusternum;i+=1)
			clusternum[unclustered_indx[Recluster_indx[i]]]=ClosestClstindx_eps[i]+1
		endfor
		
		j+=1
	while(reclusternum!=0) 
	
	extract/indx/o clusternum, unclustered_hm_indx,clusternum==-1 && weightwave >= Mthreshold 
	
	variable nunclstHM = numpnts(unclustered_hm_indx)
	
	wavestats/q clusternum
	variable ncurrentclst = V_max
	variable newclstnum = ncurrentclst+1
	
	for(j=0 ; j<nunclstHM ;j+=1)
		clusternum[unclustered_hm_indx[j]] = newclstnum
		if(newclstnum-1 >= nclusters)
			Insertpoints newclstnum-1,1,clusterseedindx
			Insertpoints/M=1 newclstnum-1,1,clusterseed
		endif
		clusterseedindx[newclstnum-1] = unclustered_hm_indx[j]
		clusterseed[][newclstnum-1] = mat[p][unclustered_hm_indx[j]]
		newclstnum+=1
	endfor
	
	ClusteronClusterNum_Wt_np(Clusternum,mat,weightwave,"_wt")
	
	wave ClusterSdv_wt,clusterupper_wt,clusterlower_wt
	killwaves recluster_indx,closestclstindx_eps, ClusterSdv_wt,clusterupper_wt,clusterlower_wt
	killwaves unclustered_indx,unclustered_hm_indx,EDindx,ClosestclstED, ClosestclstIndx
	
	string Notestr="Function used: NSSC_w1M; mat <- "+nameofwave(mat)+" ; Dmat <- "+nameofwave(Dmat)+" ; sortingwave <- "+nameofwave(sortingwave)+" ; Weightwave<- "+nameofwave(Weightwave)+" ; nclusters = "+num2istr(nclusters)+" ; eps = "+num2str(eps)
	
	Note/k clusternum "This wave stores Cluster# for each ion, with -1 indicating unclustered"
	Note clusternum Notestr
	Note/k clusterseedindx "This wave stores the index of seed for each cluster"
	Note clusterseedindx Notestr
	Note/k clusternum_org "This wave stores Cluster# for each ion before executing the cycles to recluster unclustered ions"
	Note clusternum_org Notestr
	Note/k Clusterseed "Each column of this matrix stores the thermogram of the seed used for clustering. Column 0 corresponds to Cluster1"
	Note clusterseed Notestr
	
end



//This function scans through a range of epsilon for clustering and outputs four parameters for determining optimal eps to use
//This function uses NSSC_np for clustering
Function FindOptEps_NSSC()
	
	string ymatrixstr = "cmp_matrix_sn_10a_fil"
	string weightwavestr = "cmp_mass_norm_fil"
	string wavenoisestr = "AvgNoise_ar35_fil"
	string xwavestr = "temp_c_10a"
	string suffix
	
	prompt ymatrixstr, "Select matrix of thermograms", popup WaveList("*",";","")
	prompt wavenoisestr, "Select noise wave", popup WaveList("*",";","")
	prompt weightwavestr, "Select mass contribution wave ", popup WaveList("*",";","")
	prompt xwavestr, "Select temprature wave",popup WaveList("*",";","")
	prompt suffix, "Suffix to add to the new waves"
	Doprompt "Select waves",ymatrixstr, wavenoisestr,weightwavestr,xwavestr, suffix
	
	if(V_flag)
		abort
	endif
	
	variable loweps = 1.5
	variable intervaleps = 0.1
	variable binN = 25
	variable fitendpnt = 240
	variable Mthreshold = 0.05
	
	prompt loweps, "Starting value of eps"
	prompt intervaleps, "Increment of eps in each step"
	prompt binN, "Total steps of eps"
	prompt fitendpnt, "End pnt for the down-slope fitting"
	prompt Mthreshold, "Threshold for significant ions"
	Doprompt "Enter values",loweps,intervaleps,binN,fitendpnt,Mthreshold
	
	if(V_flag)
		abort
	endif
	
	wave cmp_matrix = $ymatrixstr
	wave weightwave = $weightwavestr
	wave xwave = $xwavestr
	wave wavenoise = $wavenoisestr
	
	variable nrows=dimsize(cmp_matrix,0)
	variable ncols=dimsize(cmp_matrix,1)
	
	variable presetNc=ceil(ncols/3) //ncols/3 is used here because number larger than this will be considered as over-interpret of the data clustering
	
	make/o/d/n=(binN) $("NumofClusters_All"+suffix)  //number of clusters including one-member clusters
	make/o/d/n=(binN) $("NumofClusters_no1M"+suffix)  //number of clusters not including one-member clusters
	make/o/d/n=(binN) $("ungroupedMass"+suffix)  //unclustered mass fraction
	make/o/d/n=(binN) $("Epsilon"+suffix)
	make/o/d/n=(ncols,binN) $("ClusternumMatrix"+suffix)
	make/o/d/n=(presetNc+1,binN) $("ClustersizeMatrix"+suffix)
	make/o/d/n=(presetNc+1,binN) $("ClustermassMatrix"+suffix)	
	make/o/d/n=(presetNc,binN) $("ClusterseedIndxMatrix"+suffix)
	make/o/d/n=(binN) $("AvgInterSD"+suffix) //average distance between cluster seeds
	make/o/d/n=(binN) $("Ratio_avginterSD_eps"+suffix)
	
	make/o/d/n=(ncols)/Free tempwave
	wave numofclusterswave=$("NumofClusters_All"+suffix)
	wave numofclusters1Mwave=$("NumofClusters_no1M"+suffix)
	wave ungroupedMasswave=$("ungroupedMass"+suffix)
	wave epsilonwave=$("Epsilon"+suffix)
	wave clusternumepswave=$("ClusternumMatrix"+suffix)
	wave clustersizeepswave=$("ClustersizeMatrix"+suffix)
	wave clustermassepswave=$("ClustermassMatrix"+suffix)
	wave clusterseedindxwave=$("ClusterseedIndxMatrix"+suffix)
	wave avginterSDwave=$("AvgInterSD"+suffix)
	wave ratio_interSDwave=$("Ratio_avginterSD_eps"+suffix)

	numofclusterswave=nan
	numofclusters1Mwave=nan
	clusternumepswave=-1
	clustersizeepswave=nan
	clustermassepswave=nan
	clusterseedindxwave=nan
	avginterSDwave=nan
	ratio_interSDwave=nan
	
	variable i,j
	i=0
	variable currenteps=loweps
	variable currentClusterN
	
	CalculateED_np(cmp_matrix,"ED_fil")
	wave EuclideanD=$("ED_fil")
		
	do
		NSSC_np(cmp_matrix,EuclideanD,wavenoise,weightwave, presetNc,currenteps)
		wave clusternum
		wave clusterseedIndx
		wave clusterave_wt
		wavestats/q clusternum
		currentClusterN=V_max
		extract/Indx/o clusternum,unclustered_hm_indx,clusternum==-1 && weightwave>=Mthreshold
		variable onememberClusterN=numpnts(unclustered_hm_indx)
		variable ClstN

		ClstN = currentclusterN + onememberClusterN

		if(clstN > presetNc)
			currenteps+=intervaleps
		else
			epsilonwave[i]=currenteps
			numofclusterswave[i] = clstN
			numofclusters1Mwave[i] = clstN - onememberClusterN
			
			variable newclusternum = currentclusterN+1
			
			for(j=0;j<onememberClusterN;j+=1)
				clusternum[unclustered_hm_indx[j]]=newclusternum
				clusterseedindx[newclusternum-1] = unclustered_hm_indx[j]
				newclusternum+=1
			endfor		
			
			ClusteronClusterNum_Wt_np(Clusternum,cmp_matrix,weightwave,"_wt")
			SortClusters_2c_np(Clusterave_wt,Clusternum,clusterseedindx,xwave,"_sort")
			wave clusternum_sort
			wave clusterave_wt_sort
			wave clusterseedindx_sort
			clusternumepswave[][i]=clusternum_sort[p]			
			CountClstSize_np(Clusternum_sort,"Clustersize_sort")
			wave clustersize_sort
			CalculateClstMass_np(Clusternum_sort,weightwave,"Clustermass_sort")
			wave clustermass_sort
			ungroupedMasswave[i]=clustermass_sort[0]
			for(j=0;j<=clstN;j+=1)
				clustermassepswave[j][i]=clustermass_sort[j]
				clustersizeepswave[j][i]=clustersize_sort[j]
				if(j< clstN)
					clusterseedindxwave[j][i]=clusterseedindx_sort[j]
				elseif(j != clstN)
					clusterseedindxwave[j][i]=nan
				endif
			endfor
			
			make/o/d/n=(nrows,clstN)/Free clusterseed_matrix
			for(j=0;j<clstN;j+=1)
				clusterseed_matrix[][j]=cmp_matrix[p][clusterseedindx[j]]
			endfor
			CalculateED_np(clusterseed_matrix,"ED_clsd")
			wave waveclusterSD=$("ED_clsd")
			wavestats/q waveclusterSD
			avgintersDwave[i]=V_sum/(ClstN*(ClstN-1))
		
			i+=1
			currenteps+=intervaleps
		endif
		
		//FindTopions(Clusternum,cmp_mass,minNeighbors,"")
	while(i<binN)
	
	ratio_interSDwave=avginterSDwave/epsilonwave
	
	wave clusterseed, clusternum_org, clustersdv_wt, clusterupper_wt,clusterlower_wt
	killwaves Clusternum,Clusterseed,ClusterseedIndx,Clusternum_org,ClusterAve_wt,Clustersdv_wt,Clusterupper_wt,Clusterlower_wt,unclustered_hm_indx
	
	wave mass50loc_sort, mass75loc_sort,ED_clsd
	killwaves clusternum_sort,Clusterave_wt_sort,clusterseedindx_sort,Clustersize_sort,clustermass_sort,mass50loc_sort, mass75loc_sort,ED_clsd

	string notestr="Function FindOptEps_NSSC is used. cmp_matrix <- "+nameofwave(cmp_matrix)+" ; weightwave <- "+nameofwave(weightwave)+" ; wavenoise <- "+nameofwave(wavenoise)+" and xwave <- "+nameofwave(xwave)
	string notestr1="Parameters used: loweps = "+num2str(loweps)+" ; intervaleps = "+num2str(intervaleps)+" ; binN = "+num2str(binN)+" ; fitendpnt = "+num2str(fitendpnt)+ " ; Mthreshold = "+num2str(Mthreshold)
	
	Note/k numofclusterswave "This wave stores the number of clusters including one-member clusters for every eps scanned"
	Note numofclusterswave notestr
	Note numofclusterswave notestr1
	Note/k numofclusters1Mwave "This wave stores the number of clusters not including one-member clusters for every eps scanned"
	Note numofclusters1Mwave notestr
	Note numofclusters1Mwave notestr1
	Note/k ungroupedMasswave "This wave stores the unclustered mass fraction for every eps scanned"
	Note ungroupedMasswave notestr
	Note ungroupedMasswave notestr1
	Note/k epsilonwave "This wave stores the eps scanned"
	Note epsilonwave notestr
	Note epsilonwave notestr1
	Note/k clusternumepswave "Each column of this matrix stores the cluster# for all the ions for every eps scanned"
	Note clusternumepswave notestr
	Note clusternumepswave notestr1
	Note/k clustersizeepswave "Each column of this matrix stores the number of ions in clusters for every eps scanned"
	Note clustersizeepswave notestr
	Note clustersizeepswave notestr1
	Note/k clustermassepswave "Each column of this matrix stores the mass fraction of clusters for every eps scanned"
	Note clustermassepswave notestr
	Note clustermassepswave notestr1
	Note/k clusterseedindxwave "Each column of this matrix stores the index of seeds of the clusters for every eps scanned"
	Note clusterseedindxwave notestr
	Note clusterseedindxwave notestr1
	Note/k avginterSDwave "This wave stores the average intercluster Euclidean distance for every eps scanned"
	Note avginterSDwave notestr
	Note avginterSDwave notestr1
	Note/k ratio_interSDwave "This wave stores the ratio of average intercluster ED over epsilon for every eps scanned"
	Note ratio_interSDwave notestr
	Note ratio_interSDwave notestr1

End


Function FindOptEps_NSSC_np(cmp_matrix,loweps,intervaleps,binN,fitendpnt,Mthreshold,weightwave,xwave,wavenoise,suffix)
	wave cmp_matrix
	variable loweps,intervaleps,binN,fitendpnt,Mthreshold
	wave weightwave,xwave
	wave wavenoise
	string suffix
	
	variable nrows=dimsize(cmp_matrix,0)
	variable ncols=dimsize(cmp_matrix,1)
	
	variable presetNc=ceil(ncols/3) //ncols/3 is used here because number larger than this will be considered as over-interpret of the data clustering
	
	make/o/d/n=(binN) $("NumofClusters_All"+suffix)  //number of clusters including one-member clusters
	make/o/d/n=(binN) $("NumofClusters_no1M"+suffix)  //number of clusters not including one-member clusters
	make/o/d/n=(binN) $("ungroupedMass"+suffix)  //unclustered mass fraction
	make/o/d/n=(binN) $("Epsilon"+suffix)
	make/o/d/n=(ncols,binN) $("ClusternumMatrix"+suffix)
	make/o/d/n=(presetNc+1,binN) $("ClustersizeMatrix"+suffix)
	make/o/d/n=(presetNc+1,binN) $("ClustermassMatrix"+suffix)	
	make/o/d/n=(presetNc,binN) $("ClusterseedIndxMatrix"+suffix)
	make/o/d/n=(binN) $("AvgInterSD"+suffix) //average distance between cluster seeds
	make/o/d/n=(binN) $("Ratio_avginterSD_eps"+suffix)
	
	make/o/d/n=(ncols)/Free tempwave
	wave numofclusterswave=$("NumofClusters_All"+suffix)
	wave numofclusters1Mwave=$("NumofClusters_no1M"+suffix)
	wave ungroupedMasswave=$("ungroupedMass"+suffix)
	wave epsilonwave=$("Epsilon"+suffix)
	wave clusternumepswave=$("ClusternumMatrix"+suffix)
	wave clustersizeepswave=$("ClustersizeMatrix"+suffix)
	wave clustermassepswave=$("ClustermassMatrix"+suffix)
	wave clusterseedindxwave=$("ClusterseedIndxMatrix"+suffix)
	wave avginterSDwave=$("AvgInterSD"+suffix)
	wave ratio_interSDwave=$("Ratio_avginterSD_eps"+suffix)

	numofclusterswave=nan
	numofclusters1Mwave=nan
	clusternumepswave=-1
	clustersizeepswave=nan
	clustermassepswave=nan
	clusterseedindxwave=nan
	avginterSDwave=nan
	ratio_interSDwave=nan
	
	variable i,j
	i=0
	variable currenteps=loweps
	variable currentClusterN
	
	CalculateED_np(cmp_matrix,"ED_fil")
	wave EuclideanD=$("ED_fil")
		
	do
		NSSC_np(cmp_matrix,EuclideanD,wavenoise,weightwave, presetNc,currenteps)
		wave clusternum
		wave clusterseedIndx
		wave clusterave_wt
		wavestats/q clusternum
		currentClusterN=V_max
		extract/Indx/o clusternum,unclustered_hm_indx,clusternum==-1 && weightwave>=Mthreshold
		variable onememberClusterN=numpnts(unclustered_hm_indx)
		variable ClstN

		ClstN = currentclusterN + onememberClusterN

		if(clstN > presetNc)
			currenteps+=intervaleps
		else
			epsilonwave[i]=currenteps
			numofclusterswave[i] = clstN
			numofclusters1Mwave[i] = clstN - onememberClusterN
			
			variable newclusternum = currentclusterN+1
			
			for(j=0;j<onememberClusterN;j+=1)
				clusternum[unclustered_hm_indx[j]]=newclusternum
				clusterseedindx[newclusternum-1] = unclustered_hm_indx[j]
				newclusternum+=1
			endfor		
			
			ClusteronClusterNum_Wt_np(Clusternum,cmp_matrix,weightwave,"_wt")
			SortClusters_2c_np(Clusterave_wt,Clusternum,clusterseedindx,xwave,"_sort")
			wave clusternum_sort
			wave clusterave_wt_sort
			wave clusterseedindx_sort
			clusternumepswave[][i]=clusternum_sort[p]			
			CountClstSize_np(Clusternum_sort,"Clustersize_sort")
			wave clustersize_sort
			CalculateClstMass_np(Clusternum_sort,weightwave,"Clustermass_sort")
			wave clustermass_sort
			ungroupedMasswave[i]=clustermass_sort[0]
			for(j=0;j<=clstN;j+=1)
				clustermassepswave[j][i]=clustermass_sort[j]
				clustersizeepswave[j][i]=clustersize_sort[j]
				if(j< clstN)
					clusterseedindxwave[j][i]=clusterseedindx_sort[j]
				elseif(j != clstN)
					clusterseedindxwave[j][i]=nan
				endif
			endfor
			
			make/o/d/n=(nrows,clstN)/Free clusterseed_matrix
			for(j=0;j<clstN;j+=1)
				clusterseed_matrix[][j]=cmp_matrix[p][clusterseedindx[j]]
			endfor
			CalculateED_np(clusterseed_matrix,"ED_clsd")
			wave waveclusterSD=$("ED_clsd")
			wavestats/q waveclusterSD
			avgintersDwave[i]=V_sum/(ClstN*(ClstN-1))
		
			i+=1
			currenteps+=intervaleps
		endif
		
		//FindTopions(Clusternum,cmp_mass,minNeighbors,"")
	while(i<binN)
	
	ratio_interSDwave=avginterSDwave/epsilonwave
	
	wave clusterseed, clusternum_org, clustersdv_wt, clusterupper_wt,clusterlower_wt
	killwaves Clusternum,Clusterseed,ClusterseedIndx,Clusternum_org,ClusterAve_wt,Clustersdv_wt,Clusterupper_wt,Clusterlower_wt,unclustered_hm_indx
	
	wave mass50loc_sort, mass75loc_sort,ED_clsd
	killwaves clusternum_sort,Clusterave_wt_sort,clusterseedindx_sort,Clustersize_sort,clustermass_sort,mass50loc_sort, mass75loc_sort,ED_clsd

	string notestr="Function FindOptEps_NSSC_np is used. cmp_matrix <- "+nameofwave(cmp_matrix)+" ; weightwave <- "+nameofwave(weightwave)+" ; wavenoise <- "+nameofwave(wavenoise)+" and xwave <- "+nameofwave(xwave)
	string notestr1="Parameters used: loweps = "+num2str(loweps)+" ; intervaleps = "+num2str(intervaleps)+" ; binN = "+num2str(binN)+" ; fitendpnt = "+num2str(fitendpnt)+ " ; Mthreshold = "+num2str(Mthreshold)
	
	Note/k numofclusterswave "This wave stores the number of clusters including one-member clusters for every eps scanned"
	Note numofclusterswave notestr
	Note numofclusterswave notestr1
	Note/k numofclusters1Mwave "This wave stores the number of clusters not including one-member clusters for every eps scanned"
	Note numofclusters1Mwave notestr
	Note numofclusters1Mwave notestr1
	Note/k ungroupedMasswave "This wave stores the unclustered mass fraction for every eps scanned"
	Note ungroupedMasswave notestr
	Note ungroupedMasswave notestr1
	Note/k epsilonwave "This wave stores the eps scanned"
	Note epsilonwave notestr
	Note epsilonwave notestr1
	Note/k clusternumepswave "Each column of this matrix stores the cluster# for all the ions for every eps scanned"
	Note clusternumepswave notestr
	Note clusternumepswave notestr1
	Note/k clustersizeepswave "Each column of this matrix stores the number of ions in clusters for every eps scanned"
	Note clustersizeepswave notestr
	Note clustersizeepswave notestr1
	Note/k clustermassepswave "Each column of this matrix stores the mass fraction of clusters for every eps scanned"
	Note clustermassepswave notestr
	Note clustermassepswave notestr1
	Note/k clusterseedindxwave "Each column of this matrix stores the index of seeds of the clusters for every eps scanned"
	Note clusterseedindxwave notestr
	Note clusterseedindxwave notestr1
	Note/k avginterSDwave "This wave stores the average intercluster Euclidean distance for every eps scanned"
	Note avginterSDwave notestr
	Note avginterSDwave notestr1
	Note/k ratio_interSDwave "This wave stores the ratio of average intercluster ED over epsilon for every eps scanned"
	Note ratio_interSDwave notestr
	Note ratio_interSDwave notestr1

End


//This function calculates the unweighted average thermogram of a cluster and the upper/lower thermograms are given by Avg+/-Sdev based on clusternum
Function ClusteronClusterNum_unWt()
	
	string clusternumstr = "clusternumsort"
	string ymatrixstr = "cmp_matrix_sn_10a_fil"
	string suffix = "_unwtsort"
	
	prompt clusternumstr, "Select wave of cluster#", popup WaveList("*",";","")
	prompt ymatrixstr, "Select matrix of thermograms",popup WaveList("*",";","")
	prompt suffix, "The suffix to add to the new waves"
	Doprompt "Select waves", clusternumstr, ymatrixstr, suffix
	
	if(V_flag)
		abort
	endif
	
	wave ClusterNum = $clusternumstr
	wave ymatrix = $ymatrixstr
	
	variable nrows=dimsize(ymatrix,0)
	variable ncols=dimsize(ymatrix,1)
	wavestats/q clusternum
	variable nclusters=V_max
	 variable i,j,k
	
	if(nclusters == -1)
		nclusters = 0
	endif
	
	make/o/d/n=(nrows,nclusters) $("ClusterAve"+suffix)
	wave Clusteravewave=$("ClusterAve"+suffix)
	make/o/d/n=(nrows,nclusters) $("ClusterSdv"+suffix)
	wave ClusterSdvwave=$("ClusterSdv"+suffix)
	clusteravewave=0
	clustersdvwave=0
	variable clustersize1
	
	for(i=0;i<nclusters;i+=1)
		k=0
		extract/Free/Indx Clusternum,indexwave,clusternum==i+1
		clustersize1=numpnts(indexwave)
		
		for(j=0;j<clustersize1;j+=1)
			ClusterAvewave[][i]+=ymatrix[p][indexwave[j]]
		endfor
		ClusterAvewave[][i]/=clustersize1
		for(j=0;j<clustersize1;j+=1)
			ClusterSdvwave[][i]+=(ymatrix[p][indexwave[j]]-Clusteravewave[p][i])^2
		endfor
		if(Clustersize1!=1)
			ClusterSdvwave[][i]=sqrt(ClusterSdvwave[p][i]/(clustersize1-1))
		else
			Clustersdvwave[][i]=0
		endif
	endfor
	
	make/o/d/n=(nrows,nclusters) $("ClusterUpper"+suffix)=clusteravewave[p][q]+clustersdvwave[p][q]
	make/o/d/n=(nrows,nclusters) $("ClusterLower"+suffix)=clusteravewave[p][q]-clustersdvwave[p][q]
	
	string notestr = "Function 'ClusteronClusterNum_unWt' is used. Clusternum <- wave ' "+nameofwave(clusternum)+" ' ; ymatrix <- matrix ' "+nameofwave(ymatrix)+" ' "
	
	Note/k clusteravewave "Each column in this matrix is the unweighted average thermogram of a cluster"
	Note clusteravewave notestr
	
	Note/k clustersdvwave "Each column in this matrix is the unweighted standard deviation of a cluster"
	Note clustersdvwave notestr
	
	Note/k $("ClusterUpper"+suffix) "Each column in this matrix is the avg+1*sdev thermogram of a cluster"
	Note $("ClusterUpper"+suffix) notestr
	
	Note/k $("ClusterLower"+suffix) "Each column in this matrix is the avg-1*sdev thermogram of a cluster"
	Note $("ClusterLower"+suffix) notestr
	 
End

Function ClusteronClusterNum_unWt_np(Clusternum,ymatrix,suffix)
	
	wave ClusterNum
	wave ymatrix
	string suffix
	
	variable nrows=dimsize(ymatrix,0)
	variable ncols=dimsize(ymatrix,1)
	wavestats/q clusternum
	variable nclusters=V_max
	variable i,j,k
	 
	if(nclusters == -1)
		nclusters = 0
	endif
	
	make/o/d/n=(nrows,nclusters) $("ClusterAve"+suffix)
	wave Clusteravewave=$("ClusterAve"+suffix)
	make/o/d/n=(nrows,nclusters) $("ClusterSdv"+suffix)
	wave ClusterSdvwave=$("ClusterSdv"+suffix)
	clusteravewave=0
	clustersdvwave=0
	variable clustersize1
	
	for(i=0;i<nclusters;i+=1)
		k=0
		extract/Free/Indx Clusternum,indexwave,clusternum==i+1
		clustersize1=numpnts(indexwave)
		
		for(j=0;j<clustersize1;j+=1)
			ClusterAvewave[][i]+=ymatrix[p][indexwave[j]]
		endfor
		ClusterAvewave[][i]/=clustersize1
		for(j=0;j<clustersize1;j+=1)
			ClusterSdvwave[][i]+=(ymatrix[p][indexwave[j]]-Clusteravewave[p][i])^2
		endfor
		if(Clustersize1!=1)
			ClusterSdvwave[][i]=sqrt(ClusterSdvwave[p][i]/(clustersize1-1))
		else
			Clustersdvwave[][i]=0
		endif
	endfor
	
	make/o/d/n=(nrows,nclusters) $("ClusterUpper"+suffix)=clusteravewave[p][q]+clustersdvwave[p][q]
	make/o/d/n=(nrows,nclusters) $("ClusterLower"+suffix)=clusteravewave[p][q]-clustersdvwave[p][q]
	
	string notestr = "Function 'ClusteronClusterNum_unWt' is used. Clusternum <- wave ' "+nameofwave(clusternum)+" ' ; ymatrix <- matrix ' "+nameofwave(ymatrix)+" ' "
	
	Note/k clusteravewave "Each column in this matrix is the unweighted average thermogram of a cluster"
	Note clusteravewave notestr
	
	Note/k clustersdvwave "Each column in this matrix is the unweighted standard deviation of a cluster"
	Note clustersdvwave notestr
	
	Note/k $("ClusterUpper"+suffix) "Each column in this matrix is the avg+1*sdev thermogram of a cluster"
	Note $("ClusterUpper"+suffix) notestr
	
	Note/k $("ClusterLower"+suffix) "Each column in this matrix is the avg-1*sdev thermogram of a cluster"
	Note $("ClusterLower"+suffix) notestr
	 
End


//This function calculates the weighted average thermogram of a cluster and the upper/lower thermograms are given by Avg+/-Sdev
Function ClusteronClusterNum_Wt()

	string clusternumstr = "clusternumsort"
	string ymatrixstr = "cmp_matrix_sn_10a_fil"
	string weightwavestr = "cmp_mass_norm_fil"
	string suffix = "_wtsort"
	
	prompt clusternumstr, "Select wave of cluster#", popup WaveList("*",";","")
	prompt ymatrixstr, "Select matrix of thermograms",popup WaveList("*",";","")
	prompt weightwavestr, "Select wave for weighting", popup WaveList("*",";","")
	prompt suffix, "The suffix to add to the new waves"
	Doprompt "Select waves", clusternumstr, ymatrixstr, weightwavestr,suffix
	
	if(V_flag)
		abort
	endif
	
	wave ClusterNum = $clusternumstr
	wave ymatrix = $ymatrixstr
	wave weightwave = $weightwavestr
	
	variable nrows=dimsize(ymatrix,0)
	variable ncols=dimsize(ymatrix,1)
	wavestats/q clusternum
	variable nclusters=V_max
	 variable i,j,k
	 
	 if(nclusters == -1)
		nclusters = 0
	endif
	
	make/o/d/n=(nrows,nclusters) $("ClusterAve"+suffix)
	wave Clusteravewave=$("ClusterAve"+suffix)
	make/o/d/n=(nrows,nclusters) $("ClusterSdv"+suffix)
	wave ClusterSdvwave=$("ClusterSdv"+suffix)
	clusteravewave=0
	clustersdvwave=0
	variable clustersize1,clustersizenoZero
	variable tempsumweight
	
	for(i=0;i<nclusters;i+=1)
		k=0
		extract/Free/Indx Clusternum,indexwave,clusternum==i+1
		extract/Free weightwave,weighting,clusternum==i+1
		extract/Free weighting,weightingnoZero,weighting!=0
		clustersize1=numpnts(indexwave)
		clustersizenoZero=numpnts(weightingnoZero)
		wavestats/q weighting
		tempsumweight=V_sum
		
		for(j=0;j<clustersize1;j+=1)
			ClusterAvewave[][i]+=ymatrix[p][indexwave[j]]*weighting[j]
		endfor
		ClusterAvewave[][i]/=tempsumweight
		for(j=0;j<clustersize1;j+=1)
			ClusterSdvwave[][i]+=weighting[j]*(ymatrix[p][indexwave[j]]-Clusteravewave[p][i])^2
		endfor
		if(ClustersizenoZero==1)
			ClusterSdvwave[][i] = 0
		elseif(ClustersizenoZero==0)
			ClusterSdvwave[][i] = 0
		else
			ClusterSdvwave[][i]=sqrt(ClusterSdvwave[p][i]/((clustersizenoZero-1)/clustersizenoZero)*tempsumweight)
		endif
			
	endfor

	make/o/d/n=(nrows,nclusters) $("ClusterUpper"+suffix)=clusteravewave[p][q]+clustersdvwave[p][q]
	make/o/d/n=(nrows,nclusters) $("ClusterLower"+suffix)=clusteravewave[p][q]-clustersdvwave[p][q]	 
	 
	 string notestr = "Function 'ClusteronClusterNum_Wt' is used. Clusternum <- wave ' "+nameofwave(clusternum)+" ' ; ymatrix <- matrix ' "+nameofwave(ymatrix)+" '  and weightwave <- wave ' "+nameofwave(weightwave)+" ' "
	
	Note/k clusteravewave "Each column in this matrix is the weighted average thermogram of a cluster"
	Note clusteravewave notestr
	
	Note/k clustersdvwave "Each column in this matrix is the weighted standard deviation of a cluster"
	Note clustersdvwave notestr
	
	Note/k $("ClusterUpper"+suffix) "Each column in this matrix is the weighted avg+1*sdev thermogram of a cluster"
	Note $("ClusterUpper"+suffix) notestr
	
	Note/k $("ClusterLower"+suffix) "Each column in this matrix is the weighted avg-1*sdev thermogram of a cluster"
	Note $("ClusterLower"+suffix) notestr
	 
End

Function ClusteronClusterNum_Wt_np(Clusternum,ymatrix,weightwave,suffix)
	wave ClusterNum
	wave ymatrix
	wave weightwave
	string suffix
	
	variable nrows=dimsize(ymatrix,0)
	variable ncols=dimsize(ymatrix,1)
	wavestats/q clusternum
	variable nclusters=V_max
	 variable i,j,k
	 
	 if(nclusters == -1)
		nclusters = 0
	endif
	
	make/o/d/n=(nrows,nclusters) $("ClusterAve"+suffix)
	wave Clusteravewave=$("ClusterAve"+suffix)
	make/o/d/n=(nrows,nclusters) $("ClusterSdv"+suffix)
	wave ClusterSdvwave=$("ClusterSdv"+suffix)
	clusteravewave=0
	clustersdvwave=0
	variable clustersize1,clustersizenoZero
	variable tempsumweight
	
	for(i=0;i<nclusters;i+=1)
		k=0
		extract/Free/Indx Clusternum,indexwave,clusternum==i+1
		extract/Free weightwave,weightning,clusternum==i+1
		extract/Free weightning,weightningnoZero,weightning!=0
		clustersize1=numpnts(indexwave)
		clustersizenoZero=numpnts(weightningnoZero)
		wavestats/q weightning
		tempsumweight=V_sum
		
		for(j=0;j<clustersize1;j+=1)
			ClusterAvewave[][i]+=ymatrix[p][indexwave[j]]*weightning[j]
		endfor
		ClusterAvewave[][i]/=tempsumweight
		for(j=0;j<clustersize1;j+=1)
			ClusterSdvwave[][i]+=weightning[j]*(ymatrix[p][indexwave[j]]-Clusteravewave[p][i])^2
		endfor
		if(ClustersizenoZero==1)
			ClusterSdvwave[][i] = 0
		elseif(ClustersizenoZero==0)
			ClusterSdvwave[][i] = 0
		else
			ClusterSdvwave[][i]=sqrt(ClusterSdvwave[p][i]/((clustersizenoZero-1)/clustersizenoZero)*tempsumweight)
		endif	
	endfor

	make/o/d/n=(nrows,nclusters) $("ClusterUpper"+suffix)=clusteravewave[p][q]+clustersdvwave[p][q]
	make/o/d/n=(nrows,nclusters) $("ClusterLower"+suffix)=clusteravewave[p][q]-clustersdvwave[p][q]	 
	 
	 string notestr = "Function 'ClusteronClusterNum_Wt' is used. Clusternum <- wave ' "+nameofwave(clusternum)+" ' ; ymatrix <- matrix ' "+nameofwave(ymatrix)+" '  and weightwave <- wave ' "+nameofwave(weightwave)+" ' "
	
	Note/k clusteravewave "Each column in this matrix is the weighted average thermogram of a cluster"
	Note clusteravewave notestr
	
	Note/k clustersdvwave "Each column in this matrix is the weighted standard deviation of a cluster"
	Note clustersdvwave notestr
	
	Note/k $("ClusterUpper"+suffix) "Each column in this matrix is the weighted avg+1*sdev thermogram of a cluster"
	Note $("ClusterUpper"+suffix) notestr
	
	Note/k $("ClusterLower"+suffix) "Each column in this matrix is the weighted avg-1*sdev thermogram of a cluster"
	Note $("ClusterLower"+suffix) notestr
	 
End



//This function finds the ED between an indexed column in ymatrix and the cluster closest to it
Function FindED_ClosestCluster(indxnum,ymatrix,Avematrix)
	variable indxnum
	wave ymatrix, Avematrix
	
	variable ngroups=dimsize(Avematrix,1)
	variable nrows=dimsize(Avematrix,0)
	
	make/o/d/n=(ngroups)/Free EDwithClusters
	make/o/d/n=(nrows)/Free this,that,squared
	make/o/d/n=2 EDindx
	
	variable i
	
	for(i=0;i<ngroups;i+=1)
		this=ymatrix[p][indxnum]
		that=Avematrix[p][i]
		squared=this-that
		squared=squared^2
		wavestats/q squared
		EDwithClusters[i]=V_sum^0.5
	endfor
	wavestats/q EDwithClusters
	
	EDindx[0] = V_minloc
	EDindx[1] = V_min
	
	Note/k EDindx "This wave is a result of using function FindED_ClosestCluster()"
	Note EDindx "indxnum <- "+num2istr(indxnum)
	Note EDindx "ymatrix <- ' "+nameofwave(ymatrix)+" ' "
	Note EDindx "Avematrix <- ' "+nameofwave(Avematrix)+" ' "
	Note EDindx "The first point records the index# of the closest column in Avematrix"
	Note EDindx "The second point records the Euclidean distance to the closest column"
	
End

//This function sort the average thermograms of clusters according to  the location where 50% of mass is desorbed 
//The index of cluster seeds are sorted accordingly
Function SortClusters_1c()
	
	string Clusteravestr = "Clusterave_wt"
	string Clusternumstr = "Clusternum"
	string Clusterseedindxstr = "Clusterseedindx"
	string xwavestr = "temp_c_10a"
	string suffix = "sort"
	
	prompt clusteravestr, "Select matrix containing average thermograms", popup WaveList("*",";","")
	prompt clusternumstr, "Select wave containing cluster#", popup WaveList("*",";","")
	prompt clusterseedindxstr, "Select wave containing index of seeds", popup WaveList("*",";","")
	prompt xwavestr, "Select the temprature wave",popup WaveList("*",";","")
	prompt suffix, "The suffix to add on the new waves"
	Doprompt "Select waves", clusteravestr, clusternumstr, clusterseedindxstr, xwavestr, suffix
	
	if(V_flag)
		abort
	endif
	
	wave ClusterAve = $clusteravestr
	wave Clusternum = $clusternumstr
	wave Clusterseedindx = $clusterseedindxstr
	wave xwave = $xwavestr
	
	variable nrows=dimsize(clusterave,0)
	variable ncols=dimsize(clusterave,1)
	variable nspecies=numpnts(clusternum)
	
	make/o/d/n=(nspecies) $(nameofwave(Clusternum)+suffix)
	wave newwave=$(nameofwave(Clusternum)+suffix)
	newwave=nan
	make/o/d/n=(nrows) $(nameofwave(xwave)+"_im") //For the purpose of fitting of down slope
	wave xwave_im=$(nameofwave(xwave)+"_im")
	xwave_im=xwave
	Insertpoints nrows,1,xwave_im
	xwave_im[nrows]=xwave_im[nrows-1]+(xwave[nrows-1]-xwave[nrows-2])
	
	make/o/d/n=(ncols)/Free Mass50Loc
	Mass50Loc=nan
	make/o/d/n=(ncols)/Free GroupNum
	GroupNum=x+1
	make/o/d/n=(nrows) tempwave1
	
	variable i,j
	variable fitstartpnt
	
	for(i=0;i<ncols;i+=1)
		tempwave1=clusterave[p][i]
		wavestats/q tempwave1
		if(V_npnts==0 || V_sum==0)
			GroupNum[i]=nan
			Mass50Loc[i]=nan
		else
			Integrate tempwave1/X=xwave_im/D=temp_INT
			wave temp_INT
			for(j=0;j<nrows-1;j+=1)
				if(temp_INT[j]<=0.5*temp_INT[nrows-1] && temp_INT[j+1]>=0.5*temp_INT[nrows-1])
					Mass50Loc[i]=j
				endif
			endfor
			
		endif
	endfor
	
	extract/Free GroupNum,tempGroupnum,numtype(Groupnum)==0
	extract/Free Mass50Loc,tempMass50loc,numtype(Groupnum)==0
	sort {tempMass50Loc} tempGroupnum
	
	variable Ngroup=numpnts(tempGroupnum)
	variable count=0
	
	make/o/d/n=(Ngroup) NewGroupNum
	Newgroupnum=tempGroupNum
	
	variable Ngroup1, Ngroup2
	
	//print Ngroup
	
	for(i=0;i<nspecies;i+=1)
		if(clusternum[i]==-1)
			newwave[i]=-1
		elseif(numtype(clusternum[i])==2)
			newwave[i]=nan
		else
			for(j=0;j<Ngroup;j+=1)
				Ngroup2=j+1
				if(NewGroupNum[j]==clusternum[i])
					break
				endif
			endfor
			newwave[i]=Ngroup2
		endif
	endfor
	
	make/o/d/n=(Ngroup) $("Mass50Loc"+suffix)
	wave mass50loc_sort=$("Mass50Loc"+suffix)
	
	make/o/d/n=(nrows,Ngroup) $(nameofwave(Clusterave)+suffix)
	wave newmatrix=$(nameofwave(Clusterave)+suffix)
	make/o/d/n=(Ngroup) $(nameofwave(clusterseedindx)+suffix)
	wave newseed=$(nameofwave(clusterseedindx)+suffix)
	
	for(i=0;i<Ngroup;i+=1)
		Ngroup1=Newgroupnum[i]
		Mass50loc_sort[i]=Mass50loc[Ngroup1-1]
		Newmatrix[][i]=Clusterave[p][Ngroup1-1]
		newseed[i]=clusterseedindx[Ngroup1-1]
	endfor
	
	
	killwaves temp_INT
	
	killwaves temp_INT,mass50loc,newgroupnum, tempwave1,xwave_im
	
	string Notestr = "Function used: sortClusters_1c; clusterave <- "+nameofwave(clusterave)+" ; Clusternum <- "+nameofwave(clusternum)+" ; clusterseedindx <- "+nameofwave(clusterseedindx)+" ; xwave<- "+nameofwave(xwave)
	
	Note/k newwave "This wave stores the sorted cluster# for each ion"
	Note newwave Notestr
	Note/k Mass50Loc_sort "This wave stores the location (point) where 50% of the mass is desorbed"
	Note Mass50loc_sort Notestr
	Note/k Newmatrix "Columns of this matrix store the sorted thermograms of clusters"
	Note newmatrix Notestr
	Note/k Newseed "This wave stores the sorted index of seeds of clusters"
	Note Newseed Notestr
	
End

Function SortClusters_1c_np(Clusterave,Clusternum,Clusterseedindx,xwave,suffix)
	
	wave ClusterAve
	wave Clusternum
	wave Clusterseedindx
	wave xwave 
	string suffix
	
	string Clusteravestr = nameofwave(clusterave)
	string Clusternumstr = nameofwave(clusternum)
	string Clusterseedindxstr = nameofwave(clusterseedindx)
	string xwavestr = nameofwave(xwave)
	
	variable nrows=dimsize(clusterave,0)
	variable ncols=dimsize(clusterave,1)
	variable nspecies=numpnts(clusternum)
	
	make/o/d/n=(nspecies) $(nameofwave(Clusternum)+suffix)
	wave newwave=$(nameofwave(Clusternum)+suffix)
	newwave=nan
	make/o/d/n=(nrows) $(nameofwave(xwave)+"_im") //For the purpose of fitting of down slope
	wave xwave_im=$(nameofwave(xwave)+"_im")
	xwave_im=xwave
	Insertpoints nrows,1,xwave_im
	xwave_im[nrows]=xwave_im[nrows-1]+(xwave[nrows-1]-xwave[nrows-2])
	
	make/o/d/n=(ncols)/Free Mass50Loc
	Mass50Loc=nan
	make/o/d/n=(ncols)/Free GroupNum
	GroupNum=x+1
	make/o/d/n=(nrows) tempwave1
	
	variable i,j
	variable fitstartpnt
	
	for(i=0;i<ncols;i+=1)
		tempwave1=clusterave[p][i]
		wavestats/q tempwave1
		if(V_npnts==0 || V_sum==0)
			GroupNum[i]=nan
			Mass50Loc[i]=nan
		else
			Integrate tempwave1/X=xwave_im/D=temp_INT
			wave temp_INT
			for(j=0;j<nrows-1;j+=1)
				if(temp_INT[j]<=0.5*temp_INT[nrows-1] && temp_INT[j+1]>=0.5*temp_INT[nrows-1])
					Mass50Loc[i]=j
				endif
			endfor
			
		endif
	endfor
	
	extract/Free GroupNum,tempGroupnum,numtype(Groupnum)==0
	extract/Free Mass50Loc,tempMass50loc,numtype(Groupnum)==0
	sort {tempMass50Loc} tempGroupnum
	
	variable Ngroup=numpnts(tempGroupnum)
	variable count=0
	
	make/o/d/n=(Ngroup) NewGroupNum
	Newgroupnum=tempGroupNum
	
	variable Ngroup1, Ngroup2
	
	//print Ngroup
	
	for(i=0;i<nspecies;i+=1)
		if(clusternum[i]==-1)
			newwave[i]=-1
		elseif(numtype(clusternum[i])==2)
			newwave[i]=nan
		else
			for(j=0;j<Ngroup;j+=1)
				Ngroup2=j+1
				if(NewGroupNum[j]==clusternum[i])
					break
				endif
			endfor
			newwave[i]=Ngroup2
		endif
	endfor
	
	make/o/d/n=(Ngroup) $("Mass50Loc"+suffix)
	wave mass50loc_sort=$("Mass50Loc"+suffix)
	
	make/o/d/n=(nrows,Ngroup) $(nameofwave(Clusterave)+suffix)
	wave newmatrix=$(nameofwave(Clusterave)+suffix)
	make/o/d/n=(Ngroup) $(nameofwave(clusterseedindx)+suffix)
	wave newseed=$(nameofwave(clusterseedindx)+suffix)
	
	for(i=0;i<Ngroup;i+=1)
		Ngroup1=Newgroupnum[i]
		Mass50loc_sort[i]=Mass50loc[Ngroup1-1]
		Newmatrix[][i]=Clusterave[p][Ngroup1-1]
		newseed[i]=clusterseedindx[Ngroup1-1]
	endfor
	
	
	killwaves temp_INT
	
	killwaves temp_INT,mass50loc,newgroupnum, tempwave1,xwave_im
	
	string Notestr = "Function used: sortClusters_1c_np; clusterave <- "+nameofwave(clusterave)+" ; Clusternum <- "+nameofwave(clusternum)+" ; clusterseedindx <- "+nameofwave(clusterseedindx)+" ; xwave<- "+nameofwave(xwave)
	
	Note/k newwave "This wave stores the sorted cluster# for each ion"
	Note newwave Notestr
	Note/k Mass50Loc_sort "This wave stores the location (point) where 50% of the mass is desorbed"
	Note Mass50loc_sort Notestr
	Note/k Newmatrix "Columns of this matrix store the sorted thermograms of clusters"
	Note newmatrix Notestr
	Note/k Newseed "This wave stores the sorted index of seeds of clusters"
	Note Newseed Notestr
	
End

//This function sort the average thermograms of clusters according to 1) the location where 50% of mass is desorbed 2)the location where 75% of mass is desorbed 
//The index of cluster seeds are sorted accordingly
Function SortClusters_2c()
	
	string Clusteravestr = "Clusterave_wt"
	string Clusternumstr = "Clusternum"
	string Clusterseedindxstr = "Clusterseedindx"
	string xwavestr = "temp_c_10a"
	string suffix = "sort"
	
	prompt clusteravestr, "Select matrix containing average thermograms", popup WaveList("*",";","")
	prompt clusternumstr, "Select wave containing cluster#", popup WaveList("*",";","")
	prompt clusterseedindxstr, "Select wave containing index of seeds", popup WaveList("*",";","")
	prompt xwavestr, "Select the temprature wave",popup WaveList("*",";","")
	prompt suffix, "The suffix to add on the new waves"
	Doprompt "Select waves", clusteravestr, clusternumstr, clusterseedindxstr, xwavestr, suffix
	
	if(V_flag)
		abort
	endif
	
	wave ClusterAve = $clusteravestr
	wave Clusternum = $clusternumstr
	wave Clusterseedindx = $clusterseedindxstr
	wave xwave = $xwavestr
	
	variable nrows=dimsize(clusterave,0)
	variable ncols=dimsize(clusterave,1)
	variable nspecies=numpnts(clusternum)
	
	make/o/d/n=(nspecies) $(nameofwave(Clusternum)+suffix)
	wave newwave=$(nameofwave(Clusternum)+suffix)
	newwave=nan
	make/o/d/n=(nrows) $(nameofwave(xwave)+"_im") //For the purpose of fitting of down slope
	wave xwave_im=$(nameofwave(xwave)+"_im")
	xwave_im=xwave
	Insertpoints nrows,1,xwave_im
	xwave_im[nrows]=xwave_im[nrows-1]+(xwave[nrows-1]-xwave[nrows-2])
	
	make/o/d/n=(ncols)/Free Mass50Loc, Mass75Loc
	Mass50Loc=nan
	Mass75Loc=nan
	make/o/d/n=(ncols)/Free GroupNum
	GroupNum=x+1
	make/o/d/n=(nrows) tempwave1
	
	variable i,j
	variable fitstartpnt
	
	for(i=0;i<ncols;i+=1)
		tempwave1=clusterave[p][i]
		wavestats/q tempwave1
		if(V_npnts==0 || V_sum==0)
			GroupNum[i]=nan
			Mass50Loc[i]=nan
			Mass75Loc[i]=nan
		else
			Integrate tempwave1/X=xwave_im/D=temp_INT
			wave temp_INT
			for(j=0;j<nrows-1;j+=1)
				if(temp_INT[j]<=0.5*temp_INT[nrows-1] && temp_INT[j+1]>=0.5*temp_INT[nrows-1])
					Mass50Loc[i]=j
				endif
				if(temp_INT[j]<=0.75*temp_INT[nrows-1] && temp_INT[j+1]>=0.75*temp_INT[nrows-1])
					Mass75Loc[i]=j
				endif
			endfor
			
		endif
	endfor
	
	extract/Free GroupNum,tempGroupnum,numtype(Groupnum)==0
	extract/Free Mass50Loc,tempMass50loc,numtype(Groupnum)==0
	extract/Free Mass75Loc,tempMass75loc,numtype(Groupnum)==0
	sort {tempMass50Loc,tempMass75loc} tempGroupnum
	
	variable Ngroup=numpnts(tempGroupnum)
	variable count=0
	
	make/o/d/n=(Ngroup) NewGroupNum
	Newgroupnum=tempGroupNum
	
	variable Ngroup1, Ngroup2
	
	//print Ngroup
	
	for(i=0;i<nspecies;i+=1)
		if(clusternum[i]==-1)
			newwave[i]=-1
		elseif(numtype(clusternum[i])==2)
			newwave[i]=nan
		else
			for(j=0;j<Ngroup;j+=1)
				Ngroup2=j+1
				if(NewGroupNum[j]==clusternum[i])
					break
				endif
			endfor
			newwave[i]=Ngroup2
		endif
	endfor
	
	make/o/d/n=(Ngroup) $("Mass50Loc"+suffix)
	wave mass50loc_sort=$("Mass50Loc"+suffix)
	make/o/d/n=(Ngroup) $("Mass75Loc"+suffix)
	wave mass75loc_sort=$("Mass75Loc"+suffix)
	
	make/o/d/n=(nrows,Ngroup) $(nameofwave(Clusterave)+suffix)
	wave newmatrix=$(nameofwave(Clusterave)+suffix)
	make/o/d/n=(Ngroup) $(nameofwave(clusterseedindx)+suffix)
	wave newseed=$(nameofwave(clusterseedindx)+suffix)
	
	for(i=0;i<Ngroup;i+=1)
		Ngroup1=Newgroupnum[i]
		Mass50loc_sort[i]=Mass50loc[Ngroup1-1]
		Mass75loc_sort[i]=Mass75loc[Ngroup1-1]
		Newmatrix[][i]=Clusterave[p][Ngroup1-1]
		newseed[i]=clusterseedindx[Ngroup1-1]
	endfor
	
	
	killwaves temp_INT
	
	killwaves temp_INT,mass50loc,mass75loc, newgroupnum, tempwave1,xwave_im
	
	string Notestr = "Function used: sortClusters_2c; clusterave <- "+nameofwave(clusterave)+" ; Clusternum <- "+nameofwave(clusternum)+" ; clusterseedindx <- "+nameofwave(clusterseedindx)+" ; xwave<- "+nameofwave(xwave)
	
	Note/k newwave "This wave stores the sorted cluster# for each ion"
	Note newwave Notestr
	Note/k Mass50Loc_sort "This wave stores the location (point) where 50% of the mass is desorbed"
	Note Mass50loc_sort Notestr
	Note/k Mass75Loc_sort "This wave stores the location (point) where 75% of the mass is desorbed"
	Note Mass75loc_sort Notestr
	Note/k Newmatrix "Columns of this matrix store the sorted thermograms of clusters"
	Note newmatrix Notestr
	Note/k Newseed "This wave stores the sorted index of seeds of clusters"
	Note Newseed Notestr
	
End

Function SortClusters_2c_np(Clusterave,Clusternum,Clusterseedindx,xwave,suffix)
	
	wave ClusterAve
	wave Clusternum
	wave Clusterseedindx
	wave xwave 
	string suffix
	
	string Clusteravestr = nameofwave(clusterave)
	string Clusternumstr = nameofwave(clusternum)
	string Clusterseedindxstr = nameofwave(clusterseedindx)
	string xwavestr = nameofwave(xwave)
	
	variable nrows=dimsize(clusterave,0)
	variable ncols=dimsize(clusterave,1)
	variable nspecies=numpnts(clusternum)
	
	make/o/d/n=(nspecies) $(nameofwave(Clusternum)+suffix)
	wave newwave=$(nameofwave(Clusternum)+suffix)
	newwave=nan
	make/o/d/n=(nrows) $(nameofwave(xwave)+"_im") //For the purpose of fitting of down slope
	wave xwave_im=$(nameofwave(xwave)+"_im")
	xwave_im=xwave
	Insertpoints nrows,1,xwave_im
	xwave_im[nrows]=xwave_im[nrows-1]+(xwave[nrows-1]-xwave[nrows-2])
	
	make/o/d/n=(ncols)/Free Mass50Loc, Mass75Loc
	Mass50Loc=nan
	Mass75Loc=nan
	make/o/d/n=(ncols)/Free GroupNum
	GroupNum=x+1
	make/o/d/n=(nrows) tempwave1
	
	variable i,j
	variable fitstartpnt
	
	for(i=0;i<ncols;i+=1)
		tempwave1=clusterave[p][i]
		wavestats/q tempwave1
		if(V_npnts==0 || V_sum==0)
			GroupNum[i]=nan
			Mass50Loc[i]=nan
			Mass75Loc[i]=nan
		else
			Integrate tempwave1/X=xwave_im/D=temp_INT
			wave temp_INT
			for(j=0;j<nrows-1;j+=1)
				if(temp_INT[j]<=0.5*temp_INT[nrows-1] && temp_INT[j+1]>=0.5*temp_INT[nrows-1])
					Mass50Loc[i]=j
				endif
				if(temp_INT[j]<=0.75*temp_INT[nrows-1] && temp_INT[j+1]>=0.75*temp_INT[nrows-1])
					Mass75Loc[i]=j
				endif
			endfor
			
		endif
	endfor
	
	extract/Free GroupNum,tempGroupnum,numtype(Groupnum)==0
	extract/Free Mass50Loc,tempMass50loc,numtype(Groupnum)==0
	extract/Free Mass75Loc,tempMass75loc,numtype(Groupnum)==0
	sort {tempMass50Loc,tempMass75loc} tempGroupnum
	
	variable Ngroup=numpnts(tempGroupnum)
	variable count=0
	
	make/o/d/n=(Ngroup) NewGroupNum
	Newgroupnum=tempGroupNum
	
	variable Ngroup1, Ngroup2
	
	//print Ngroup
	
	for(i=0;i<nspecies;i+=1)
		if(clusternum[i]==-1)
			newwave[i]=-1
		elseif(numtype(clusternum[i])==2)
			newwave[i]=nan
		else
			for(j=0;j<Ngroup;j+=1)
				Ngroup2=j+1
				if(NewGroupNum[j]==clusternum[i])
					break
				endif
			endfor
			newwave[i]=Ngroup2
		endif
	endfor
	
	make/o/d/n=(Ngroup) $("Mass50Loc"+suffix)
	wave mass50loc_sort=$("Mass50Loc"+suffix)
	make/o/d/n=(Ngroup) $("Mass75Loc"+suffix)
	wave mass75loc_sort=$("Mass75Loc"+suffix)
	
	make/o/d/n=(nrows,Ngroup) $(nameofwave(Clusterave)+suffix)
	wave newmatrix=$(nameofwave(Clusterave)+suffix)
	make/o/d/n=(Ngroup) $(nameofwave(clusterseedindx)+suffix)
	wave newseed=$(nameofwave(clusterseedindx)+suffix)
	
	for(i=0;i<Ngroup;i+=1)
		Ngroup1=Newgroupnum[i]
		Mass50loc_sort[i]=Mass50loc[Ngroup1-1]
		Mass75loc_sort[i]=Mass75loc[Ngroup1-1]
		Newmatrix[][i]=Clusterave[p][Ngroup1-1]
		newseed[i]=clusterseedindx[Ngroup1-1]
	endfor
	
	
	killwaves temp_INT
	
	killwaves temp_INT,mass50loc,mass75loc, newgroupnum, tempwave1,xwave_im
	
	string Notestr = "Function used: sortClusters_2c; clusterave <- "+nameofwave(clusterave)+" ; Clusternum <- "+nameofwave(clusternum)+" ; clusterseedindx <- "+nameofwave(clusterseedindx)+" ; xwave<- "+nameofwave(xwave)
	
	Note/k newwave "This wave stores the sorted cluster# for each ion"
	Note newwave Notestr
	Note/k Mass50Loc_sort "This wave stores the location (point) where 50% of the mass is desorbed"
	Note Mass50loc_sort Notestr
	Note/k Mass75Loc_sort "This wave stores the location (point) where 75% of the mass is desorbed"
	Note Mass75loc_sort Notestr
	Note/k Newmatrix "Columns of this matrix store the sorted thermograms of clusters"
	Note newmatrix Notestr
	Note/k Newseed "This wave stores the sorted index of seeds of clusters"
	Note Newseed Notestr
	
End


//This function sort the average thermograms of clusters according to 1) the location where 50% of mass is desorbed 2) down slope
//The index of cluster seeds are sorted accordingly
Function SortClusters_2c_slope()
	
	string Clusteravestr = "Clusterave_wt"
	string Clusternumstr = "Clusternum"
	string Clusterseedindxstr = "Clusterseedindx"
	string xwavestr = "temp_c_10a"
	variable fitendpnt = 240
	string suffix = "sort"
	
	prompt clusteravestr, "Select matrix containing average thermograms", popup WaveList("*",";","")
	prompt clusternumstr, "Select wave containing cluster#", popup WaveList("*",";","")
	prompt clusterseedindxstr, "Select wave containing index of seeds", popup WaveList("*",";","")
	prompt xwavestr, "Select the temprature wave",popup WaveList("*",";","")
	prompt suffix, "The suffix to add on the new waves"
	prompt fitendpnt, "The ending point of down slope fitting"
	Doprompt "Select waves", clusteravestr, clusternumstr, clusterseedindxstr, xwavestr, suffix, fitendpnt
	
	if(V_flag)
		abort
	endif
	
	wave ClusterAve = $clusteravestr
	wave Clusternum = $clusternumstr
	wave Clusterseedindx = $clusterseedindxstr
	wave xwave = $xwavestr
	
	variable nrows=dimsize(clusterave,0)
	variable ncols=dimsize(clusterave,1)
	variable nspecies=numpnts(clusternum)
	
	make/o/d/n=(nspecies) $(nameofwave(Clusternum)+suffix)
	wave newwave=$(nameofwave(Clusternum)+suffix)
	newwave=nan
	make/o/d/n=(nrows) $(nameofwave(xwave)+"_im") //For the purpose of fitting of down slope
	wave xwave_im=$(nameofwave(xwave)+"_im")
	xwave_im=xwave
	Insertpoints nrows,1,xwave_im
	xwave_im[nrows]=xwave_im[nrows-1]+(xwave[nrows-1]-xwave[nrows-2])
	
	make/o/d/n=(ncols)/Free DecaySlope
	make/o/d/n=(ncols)/Free Mass50Loc
	DecaySlope=nan
	Mass50Loc=nan
	make/o/d/n=(ncols)/Free GroupNum
	GroupNum=x+1
	make/o/d/n=(ncols)/Free DecaySlope_inv
	DecaySlope_inv=nan
	make/o/d/n=(ncols,3)/Free ExpfitCoef,ExpfitSigma
	ExpfitCoef=nan
	ExpfitSigma=nan
	make/o/d/n=(nrows) tempwave1
	
	variable i,j
	variable fitstartpnt
	
	for(i=0;i<ncols;i+=1)
		tempwave1=clusterave[p][i]
		MovingAverageWave_np(tempwave1,35)
		wave tempwave1_smt35
		wavestats/q tempwave1
		if(V_npnts==0 || V_sum==0)
			DecaySlope[i]=nan
			GroupNum[i]=nan
			Mass50Loc[i]=nan
		else
			wavestats/q tempwave1_smt35
			fitstartpnt = V_maxloc+17
			Integrate tempwave1/X=xwave_im/D=temp_INT
			wave temp_INT
			for(j=0;j<nrows-1;j+=1)
				if(temp_INT[j]<=0.5*temp_INT[nrows-1] && temp_INT[j+1]>=0.5*temp_INT[nrows-1])
					Mass50Loc[i]=j
				endif
			endfor
			if(fitstartpnt < fitendpnt-53)
				CurveFit/Q/M=2/W=0 exp, tempwave1[fitstartpnt,fitendpnt]/X=xwave[fitstartpnt,fitendpnt]/D
				wave W_coef,W_sigma
				ExpfitCoef[i][]=W_coef[q]
				ExpfitSigma[i][]=W_sigma[q]
			elseif(fitstartpnt < fitendpnt -3)
				
			else
				ExpfitCoef[i][]=nan
				ExpfitSigma[i][]=nan
			endif
			DecaySlope[i]=ExpfitCoef[i][2]
			
		endif
	endfor
	DecaySlope_inv=1/DecaySlope
	
	extract/Free GroupNum,tempGroupnum,numtype(Groupnum)==0
	extract/Free Mass50Loc,tempMass50loc,numtype(Groupnum)==0
	extract/Free DecaySlope_inv,tempDecaySlope,numtype(Groupnum)==0
	sort {tempMass50Loc,tempDecaySlope} tempGroupnum
	
	variable Ngroup=numpnts(tempGroupnum)
	variable count=0
	
	make/o/d/n=(Ngroup) NewGroupNum
	Newgroupnum=tempGroupNum
	
	variable Ngroup1, Ngroup2
	
	//print Ngroup
	
	for(i=0;i<nspecies;i+=1)
		if(clusternum[i]==-1)
			newwave[i]=-1
		elseif(numtype(clusternum[i])==2)
			newwave[i]=nan
		else
			for(j=0;j<Ngroup;j+=1)
				Ngroup2=j+1
				if(NewGroupNum[j]==clusternum[i])
					break
				endif
			endfor
			newwave[i]=Ngroup2
		endif
	endfor
	
	make/o/d/n=(Ngroup) $("Mass50Loc"+suffix)
	wave mass50loc_sort=$("Mass50Loc"+suffix)
	
	make/o/d/n=(nrows,Ngroup) $(nameofwave(Clusterave)+suffix)
	wave newmatrix=$(nameofwave(Clusterave)+suffix)
	make/o/d/n=(Ngroup) $(nameofwave(clusterseedindx)+suffix)
	wave newseed=$(nameofwave(clusterseedindx)+suffix)
	
	for(i=0;i<Ngroup;i+=1)
		Ngroup1=Newgroupnum[i]
		Mass50loc_sort[i]=Mass50loc[Ngroup1-1]
		Newmatrix[][i]=Clusterave[p][Ngroup1-1]
		newseed[i]=clusterseedindx[Ngroup1-1]
	endfor
	
	
	killwaves temp_INT
	
	wave W_sigma,M_Covar,W_coef,fit_tempwave1
	killwaves temp_INT,decayslope,mass50loc,expfitcoef,expfitsigma,newgroupnum,W_coef,M_covar,W_sigma,fit_tempwave1, tempwave1, tempwave1_smt35,xwave_im
	
	string Notestr = "Function used: sortClusters_2c_slope; clusterave <- "+nameofwave(clusterave)+" ; Clusternum <- "+nameofwave(clusternum)+" ; clusterseedindx <- "+nameofwave(clusterseedindx)+" ; xwave<- "+nameofwave(xwave)+" ; fitendpnt = "+num2istr(fitendpnt)
	
	Note/k newwave "This wave stores the sorted cluster# for each ion"
	Note newwave Notestr
	Note/k Mass50Loc_sort "This wave stores the location (point) where 50% of the mass is desorbed"
	Note Mass50loc_sort Notestr
	Note/k Newmatrix "Columns of this matrix store the sorted thermograms of clusters"
	Note newmatrix Notestr
	Note/k Newseed "This wave stores the sorted index of seeds of clusters"
	Note Newseed Notestr
	
End


Function SortClusters_2c_slope_np(Clusterave,Clusternum,Clusterseedindx,xwave,fitendpnt,suffix)
	
	wave ClusterAve, Clusternum, Clusterseedindx, xwave
	variable fitendpnt
	string suffix
	
	variable nrows=dimsize(clusterave,0)
	variable ncols=dimsize(clusterave,1)
	variable nspecies=numpnts(clusternum)
	
	make/o/d/n=(nspecies) $(nameofwave(Clusternum)+suffix)
	wave newwave=$(nameofwave(Clusternum)+suffix)
	newwave=nan
	make/o/d/n=(nrows) $(nameofwave(xwave)+"_im") //For the purpose of fitting of down slope
	wave xwave_im=$(nameofwave(xwave)+"_im")
	xwave_im=xwave
	Insertpoints nrows,1,xwave_im
	xwave_im[nrows]=xwave_im[nrows-1]+(xwave[nrows-1]-xwave[nrows-2])
	
	make/o/d/n=(ncols)/Free DecaySlope
	make/o/d/n=(ncols)/Free Mass50Loc
	DecaySlope=nan
	Mass50Loc=nan
	make/o/d/n=(ncols)/Free GroupNum
	GroupNum=x+1
	make/o/d/n=(ncols)/Free DecaySlope_inv
	DecaySlope_inv=nan
	make/o/d/n=(ncols,3)/Free ExpfitCoef,ExpfitSigma
	ExpfitCoef=nan
	ExpfitSigma=nan
	make/o/d/n=(nrows) tempwave1
	
	variable i,j
	variable fitstartpnt
	
	for(i=0;i<ncols;i+=1)
		tempwave1=clusterave[p][i]
		MovingAverageWave_np(tempwave1,35)
		wave tempwave1_smt35
		wavestats/q tempwave1
		if(V_npnts==0 || V_sum==0)
			DecaySlope[i]=nan
			GroupNum[i]=nan
			Mass50Loc[i]=nan
		else
			wavestats/q tempwave1_smt35
			fitstartpnt = V_maxloc+17
			Integrate tempwave1/X=xwave_im/D=temp_INT
			wave temp_INT
			for(j=0;j<nrows-1;j+=1)
				if(temp_INT[j]<=0.5*temp_INT[nrows-1] && temp_INT[j+1]>=0.5*temp_INT[nrows-1])
					Mass50Loc[i]=j
				endif
			endfor
			if(fitstartpnt < fitendpnt-3)
				CurveFit/Q/M=2/W=2 exp, tempwave1[fitstartpnt,fitendpnt]/X=xwave[fitstartpnt,fitendpnt]/D
				wave W_coef,W_sigma
				ExpfitCoef[i][]=W_coef[q]
				ExpfitSigma[i][]=W_sigma[q]
			else
				ExpfitCoef[i][]=nan
				ExpfitSigma[i][]=nan
			endif
			DecaySlope[i]=ExpfitCoef[i][2]
			
		endif
	endfor
	DecaySlope_inv=1/DecaySlope
	
	extract/Free GroupNum,tempGroupnum,numtype(Groupnum)==0
	extract/Free Mass50Loc,tempMass50loc,numtype(Groupnum)==0
	extract/Free DecaySlope_inv,tempDecaySlope,numtype(Groupnum)==0
	sort {tempMass50Loc,tempDecaySlope} tempGroupnum
	
	variable Ngroup=numpnts(tempGroupnum)
	variable count=0
	
	make/o/d/n=(Ngroup) NewGroupNum
	Newgroupnum=tempGroupNum
	
	variable Ngroup1, Ngroup2
	
	//print Ngroup
	
	for(i=0;i<nspecies;i+=1)
		if(clusternum[i]==-1)
			newwave[i]=-1
		elseif(numtype(clusternum[i])==2)
			newwave[i]=nan
		else
			for(j=0;j<Ngroup;j+=1)
				Ngroup2=j+1
				if(NewGroupNum[j]==clusternum[i])
					break
				endif
			endfor
			newwave[i]=Ngroup2
		endif
	endfor
	
	make/o/d/n=(Ngroup) $("Mass50Loc"+suffix)
	wave mass50loc_sort=$("Mass50Loc"+suffix)
	
	make/o/d/n=(nrows,Ngroup) $(nameofwave(Clusterave)+suffix)
	wave newmatrix=$(nameofwave(Clusterave)+suffix)
	make/o/d/n=(Ngroup) $(nameofwave(clusterseedindx)+suffix)
	wave newseed=$(nameofwave(clusterseedindx)+suffix)
	
	for(i=0;i<Ngroup;i+=1)
		Ngroup1=Newgroupnum[i]
		Mass50loc_sort[i]=Mass50loc[Ngroup1-1]
		Newmatrix[][i]=Clusterave[p][Ngroup1-1]
		newseed[i]=clusterseedindx[Ngroup1-1]
	endfor
	
	
	killwaves temp_INT
	
	wave W_sigma,M_Covar,W_coef,fit_tempwave1
	killwaves temp_INT,decayslope,mass50loc,expfitcoef,expfitsigma,newgroupnum,W_coef,M_covar,W_sigma,fit_tempwave1, tempwave1, tempwave1_smt35,xwave_im
	
	string Notestr = "Function used: sortClusters_2c_slope_np; clusterave <- "+nameofwave(clusterave)+" ; Clusternum <- "+nameofwave(clusternum)+" ; clusterseedindx <- "+nameofwave(clusterseedindx)+" ; xwave<- "+nameofwave(xwave)+" ; fitendpnt = "+num2istr(fitendpnt)
	
	Note/k newwave "This wave stores the sorted cluster# for each ion"
	Note newwave Notestr
	Note/k Mass50Loc_sort "This wave stores the location (point) where 50% of the mass is desorbed"
	Note Mass50loc_sort Notestr
	Note/k Newmatrix "Columns of this matrix store the sorted thermograms of clusters"
	Note newmatrix Notestr
	Note/k Newseed "This wave stores the sorted index of seeds of clusters"
	Note Newseed Notestr
	
End

//This function count the number of ions in each cluster, including the number of unclustered ions
Function CountClstSize()

	string clusternumstr = "Clusternumsort"
	string newwavestr = "Clustersizesort"
	
	prompt clusternumstr, "Select wave of cluster#", popup WaveList("*",";","")
	prompt newwavestr, "Name of the new wave"
	Doprompt "Select waves",clusternumstr,newwavestr
	
	if(V_flag)
		abort
	endif
	
	wave clusternum = $clusternumstr
	
	wavestats/q clusternum
	variable Groupnum=V_max
	variable nspecies=numpnts(clusternum)
	make/o/d/n=(Groupnum+1) $newwavestr
	wave newwave=$newwavestr
	
	variable i,j,count
	count=0
	
	for(i=0;i<groupnum;i+=1)
		for(j=0;j<nspecies;j+=1)
			if(clusternum[j]==i+1)
				count+=1
			endif
		endfor
		newwave[i+1]=count
		count=0
	endfor
	
	for(j=0;j<nspecies;j+=1)
		if(clusternum[j]==-1)
				count+=1
		endif
	endfor
	newwave[0]=count
	
	Note/k newwave "This wave stores the number of ions of clusters, point 0 -> unclustered, point 1-x -> Cluster 1-x"
	Note newwave "Function CountClstSize_np is used. Clusternum <- "+nameofwave(clusternum)
	
End


Function CountClstSize_np(Clusternum,newwavestr)
	wave clusternum
	string newwavestr
	
	wavestats/q clusternum
	variable Groupnum=V_max
	variable nspecies=numpnts(clusternum)
	make/o/d/n=(Groupnum+1) $newwavestr
	wave newwave=$newwavestr
	
	variable i,j,count
	count=0
	
	for(i=0;i<groupnum;i+=1)
		for(j=0;j<nspecies;j+=1)
			if(clusternum[j]==i+1)
				count+=1
			endif
		endfor
		newwave[i+1]=count
		count=0
	endfor
	
	for(j=0;j<nspecies;j+=1)
		if(clusternum[j]==-1)
				count+=1
		endif
	endfor
	newwave[0]=count
	
	Note/k newwave "This wave stores the number of ions of clusters, point 0 -> unclustered, point 1-x -> Cluster 1-x"
	Note newwave "Function CountClstSize_np is used. Clusternum <- "+nameofwave(clusternum)
	
End

//This function calculates the mass fractional contribution of each cluster by adding up the mass fra of ions in each cluster
Function CalculateClstMass()

	string clusternumstr = "Clusternumsort"
	string cmpmassstr = "cmp_mass_norm_fil"
	string newwavestr = "Clustermasssort"
	
	prompt clusternumstr, "Select wave of cluster#", popup WaveList("*",";","")
	prompt cmpmassstr, "Select wave of mass contribution", popup WaveList("*",";","")
	prompt newwavestr, "Name of the new wave"
	Doprompt "Select waves",clusternumstr,cmpmassstr, newwavestr
	
	if(V_flag)
		abort
	endif
	
	wave clusternum = $clusternumstr
	wave Cmpmass = $cmpmassstr
	
	variable nrows=numpnts(clusternum)
	variable i,j,count
	wavestats/q clusternum
	variable nclusters=V_max+1
	
	make/o/d/n=(nclusters) $newwavestr
	wave newwave=$newwavestr
	
	for(i=0;i<nclusters;i+=1)
		count=0
		if(i==0)
			for(j=0;j<nrows;j+=1)
				if(clusternum[j]==-1)
					count+=Cmpmass[j]
				endif
			endfor
		else
			for(j=0;j<nrows;j+=1)
				if(clusternum[j]==i)
					count+=Cmpmass[j]
				endif
			endfor
		endif
		newwave[i]=count
	endfor
	
	Note/k newwave "This wave stores the mass fractional contribution of clusters, point 0 -> unclustered, point 1-x -> Cluster 1-x"
	Note newwave "Function CalculateClstMass_np is used. Clusternum <- "+nameofwave(clusternum)+" and Cmpmass <- "+nameofwave(CmpMass)
	 
End

Function CalculateClstMass_np(Clusternum,CmpMass,newwavestr)
	wave clusternum,Cmpmass
	string newwavestr
	
	variable nrows=numpnts(clusternum)
	variable i,j,count
	wavestats/q clusternum
	variable nclusters=V_max+1
	
	make/o/d/n=(nclusters) $newwavestr
	wave newwave=$newwavestr
	
	for(i=0;i<nclusters;i+=1)
		count=0
		if(i==0)
			for(j=0;j<nrows;j+=1)
				if(clusternum[j]==-1)
					count+=Cmpmass[j]
				endif
			endfor
		else
			for(j=0;j<nrows;j+=1)
				if(clusternum[j]==i)
					count+=Cmpmass[j]
				endif
			endfor
		endif
		newwave[i]=count
	endfor
	
	Note/k newwave "This wave stores the mass fractional contribution of clusters, point 0 -> unclustered, point 1-x -> Cluster 1-x"
	Note newwave "Function CalculateClstMass_np is used. Clusternum <- "+nameofwave(clusternum)+" and Cmpmass <- "+nameofwave(CmpMass)
	 
End


//This function calculates the average properties of each cluster, such as average H/C, O/C
Function AverageClstProperties(weighted)
	variable weighted // 0 = no, 1 = mass weighted

	string clusternumstr = "Clusternumsort"
	string clustersizestr = "Clustersizesort"
	string weight_str = "no"

	prompt clusternumstr, "Select wave of cluster#", popup WaveList("*",";","")
	prompt clustersizestr, "Select wave of cluster size", popup WaveList("*",";","")
	prompt weight_str, "Weight by mass contributions?", popup "yes;no"
	Doprompt "Select waves",clusternumstr,clustersizestr, weight_Str
	
	if(V_flag)
		abort
	endif

	wave clusternumwave = $clusternumstr //-1 for unclustered ions and 1-x as cluster# for clustered ions
	wave clustersizewave = $clustersizestr //By default, the clustersize wave contains the number of unclustered ions in point 0

	//Make sure that the following waves already exists in the current folder
	
	wave Molarmass_cmp_fil
	wave N_num_fil
	wave C_num_fil
	wave H_num_fil
	wave O_num_fil
	wave OSc_fil
	wave OC_fil
	wave HC_fil
	wave Cmp_mass_fil
	
		
	if(waveexists(Molarmass_cmp_fil)==0 || waveexists(N_num_fil) ==0 || waveexists(C_num_fil) == 0|| waveexists(H_num_fil) == 0|| waveexists(O_num_fil) == 0|| waveexists(HC_fil) == 0|| waveexists(OC_fil) == 0|| waveexists(OSc_fil) == 0|| waveexists(Cmp_mass_fil) == 0)
		DoAlert 0, "Please make sure that waves Molarmass_cmp_fil, N_num_fil, C_num_fil, H_num_fil, O_num_fil, OSc_fil, OC_fil, HC_fil, Cmp_mass_fil are in the current folder"
		abort
	endif
	
	variable nrows = numpnts(clusternumwave)
	
	wavestats/q clusternumwave
	variable nclusters = V_max
	
	if(nclusters!=numpnts(clustersizewave)-1)
		print "Please make sure the clusternum wave and clustersize wave are results from the same clustering"
		abort
	endif
	
	variable i,j,k,clustercount
	variable numcount,tempsum
	
	make/o/d/n=(nclusters,2) Molarmass_clst_avg,OSc_clst_avg,OC_clst_avg,HC_clst_avg,C_num_clst_avg,N_num_clst_avg,O_num_clst_avg,H_num_clst_avg,Mass_clst_avg
	make/o/d/n=(nclusters) Mass_clst_sum
	make/o/t/n=(nclusters) MolecularFormula_clst_avg // a text wave with the average molecular formula for each cluster
	
	string MFstr = ""
	variable atom_ceil, atom_floor, atom_dec // for determining average molecular formula
//	variable weighted // 0 = no, 1 = weight by mass
	if(stringmatch(weight_str,"no")==1)
		weighted = 0
	else
		weighted = 1
	endif
	variable V_avg_weighted
	
	for(i=0;i<nclusters;i+=1)
		MFstr = ""
		make/o/d/n=(clustersizewave[i+1])/Free tempwave,tempwave1,tempwave2,tempwave3,tempwave4,tempwave5,tempwave6,tempwave7,tempwave8
		numcount=0
		if(weighted==0) // no weighting
			for(j=0;j<nrows;j+=1)
				if(clusternumwave[j]==i+1)
					tempwave[numcount]=Molarmass_cmp_fil[j]
					tempwave1[numcount]=OSc_fil[j]
					tempwave2[numcount]=OC_fil[j]
					tempwave3[numcount]=HC_fil[j]
					tempwave4[numcount]=C_num_fil[j]
					tempwave5[numcount]=N_num_fil[j]
					tempwave6[numcount]=Cmp_mass_fil[j]
					tempwave7[numcount]=O_num_fil[j]
					tempwave8[numcount]=H_num_fil[j]
					numcount+=1
				endif
			endfor
		
			// mass
			wavestats/q tempwave6
			mass_clst_avg[i][0]=V_avg
			mass_clst_avg[i][1]=V_sdev
			Mass_clst_sum[i]=V_sum
			wavestats/q tempwave
			Molarmass_clst_avg[i][0]=V_avg
			Molarmass_clst_avg[i][1]=V_sdev
			wavestats/q tempwave1
			OSc_clst_avg[i][0]=V_avg
			OSc_clst_avg[i][1]=V_sdev
			wavestats/q tempwave2
			OC_clst_avg[i][0]=V_avg
			OC_clst_avg[i][1]=V_sdev
			wavestats/q tempwave3
			HC_clst_avg[i][0]=V_avg
			HC_clst_avg[i][1]=V_sdev
			// Carbon atoms
			wavestats/q tempwave4
			C_num_clst_avg[i][0]=V_avg
			C_num_clst_avg[i][1]=V_sdev
			atom_ceil = ceil(V_avg)
			atom_floor = floor(V_avg)
			atom_dec = floor(10*(atom_ceil - V_avg))
			MFstr += "C" + num2istr(atom_floor) + "." + num2istr(atom_dec)
			// H atoms
			wavestats/q tempwave8
			H_num_clst_avg[i][0]=V_avg
			H_num_clst_avg[i][1]=V_sdev
			atom_ceil = ceil(V_avg)
			atom_floor = floor(V_avg)
			atom_dec = floor(10*(atom_ceil - V_avg))
			MFstr += "H" + num2istr(atom_floor) + "." + num2istr(atom_dec)
			// O atoms
			wavestats/q tempwave7
			O_num_clst_avg[i][0]=V_avg
			O_num_clst_avg[i][1]=V_sdev
			atom_ceil = ceil(V_avg)
			atom_floor = floor(V_avg)
			atom_dec = floor(10*(atom_ceil - V_avg))
			MFstr += "O" + num2istr(atom_floor) + "." + num2istr(atom_dec)
			// N atoms
			wavestats/q tempwave5
			N_num_clst_avg[i][0]=V_avg
			N_num_clst_avg[i][1]=V_sdev
			atom_ceil = ceil(V_avg)
			atom_floor = floor(V_avg)
			atom_dec = floor(10*(atom_ceil - V_avg))
			MFstr += "N" + num2istr(atom_floor) + "." + num2istr(atom_dec)
			
			MolecularFormula_clst_avg[i] = MFstr
		else // weighted by mass
			for(j=0;j<nrows;j+=1)
				if(clusternumwave[j]==i+1)
					tempwave[numcount]=Molarmass_cmp_fil[j] * Cmp_mass_fil[j]
					tempwave1[numcount]=OSc_fil[j]
					tempwave2[numcount]=OC_fil[j]
					tempwave3[numcount]=HC_fil[j]
					tempwave4[numcount]=C_num_fil[j] * Cmp_mass_fil[j]
					tempwave5[numcount]=N_num_fil[j] * Cmp_mass_fil[j]
					tempwave6[numcount]=Cmp_mass_fil[j]
					tempwave7[numcount]=O_num_fil[j] * Cmp_mass_fil[j]
					tempwave8[numcount]=H_num_fil[j] * Cmp_mass_fil[j]
					numcount+=1
				endif
			endfor
		
			// mass
			wavestats/q tempwave6
			mass_clst_avg[i][0]=V_avg
			mass_clst_avg[i][1]=V_sdev
			Mass_clst_sum[i]=V_sum
			wavestats/q tempwave
			Molarmass_clst_avg[i][0]=V_avg/mass_clst_avg[i][0]
			Molarmass_clst_avg[i][1]=V_sdev/mass_clst_avg[i][0]
			wavestats/q tempwave1
			OSc_clst_avg[i][0]=V_avg
			OSc_clst_avg[i][1]=V_sdev
			wavestats/q tempwave2
			OC_clst_avg[i][0]=V_avg
			OC_clst_avg[i][1]=V_sdev
			wavestats/q tempwave3
			HC_clst_avg[i][0]=V_avg
			HC_clst_avg[i][1]=V_sdev
			// Carbon atoms
			wavestats/q tempwave4
			C_num_clst_avg[i][0]=V_avg/mass_clst_avg[i][0]
			V_avg_weighted = V_avg/mass_clst_avg[i][0]
			C_num_clst_avg[i][1]=V_sdev/mass_clst_avg[i][0]
			atom_ceil = ceil(V_avg_weighted)
			atom_floor = floor(V_avg_weighted)
			atom_dec = round(10*V_avg_weighted)-10*floor(V_avg_weighted) // floor(10*(1-(atom_ceil - V_avg_weighted)))
			if(atom_dec==10)
				atom_floor +=1 
				atom_dec = 0
			endif
			MFstr += "C" + num2istr(atom_floor) + "." + num2istr(atom_dec)
			// H atoms
			wavestats/q tempwave8
			H_num_clst_avg[i][0]=V_avg/mass_clst_avg[i][0]
			V_avg_weighted = V_avg/mass_clst_avg[i][0]
			H_num_clst_avg[i][1]=V_sdev/mass_clst_avg[i][0]
			atom_ceil = ceil(V_avg_weighted)
			atom_floor = floor(V_avg_weighted)
			atom_dec = round(10*V_avg_weighted)-10*floor(V_avg_weighted) // floor(10*(1-(atom_ceil - V_avg_weighted)))
			if(atom_dec==10)
				atom_floor +=1 
				atom_dec = 0
			endif
			MFstr += "H" + num2istr(atom_floor) + "." + num2istr(atom_dec)
			// O atoms
			wavestats/q tempwave7
			O_num_clst_avg[i][0]=V_avg/mass_clst_avg[i][0]
			V_avg_weighted = V_avg/mass_clst_avg[i][0]
			O_num_clst_avg[i][1]=V_sdev/mass_clst_avg[i][0]
			atom_ceil = ceil(V_avg_weighted)
			atom_floor = floor(V_avg_weighted)
			atom_dec = round(10*V_avg_weighted)-10*floor(V_avg_weighted) // floor(10*(1-(atom_ceil - V_avg_weighted)))
			if(atom_dec==10)
				atom_floor +=1 
				atom_dec = 0
			endif
			MFstr += "O" + num2istr(atom_floor) + "." + num2istr(atom_dec)
			// N atoms
			wavestats/q tempwave5
			N_num_clst_avg[i][0]=V_avg/mass_clst_avg[i][0]
			V_avg_weighted = V_avg/mass_clst_avg[i][0]
			N_num_clst_avg[i][1]=V_sdev/mass_clst_avg[i][0]
			atom_ceil = ceil(V_avg_weighted)
			atom_floor = floor(V_avg_weighted)
			atom_dec = round(10*V_avg_weighted)-10*floor(V_avg_weighted) // floor(10*(1-(atom_ceil - V_avg_weighted)))
			if(atom_dec==10)
				atom_floor +=1 
				atom_dec = 0
			endif
			MFstr += "N" + num2istr(atom_floor) + "." + num2istr(atom_dec)
			
			MolecularFormula_clst_avg[i] = MFstr
		endif
		// Get O:C, H:C, etc. 
		OC_clst_avg[][0] = O_num_clst_avg[p][0]/C_num_clst_avg[p][0]
		HC_clst_avg[][0] = H_num_clst_avg[p][0]/C_num_clst_avg[p][0]
	endfor
	
	Note/k Molarmass_clst_avg "This wave stores the average molarmass of clusters (column0) and sdev (column1) according to wave "+nameofwave(clusternumwave)
	Note/k OSc_clst_avg "This wave stores the average OSc of clusters (column0) and sdev (column1) according to wave "+nameofwave(clusternumwave)
	Note/k OC_clst_avg "This wave stores the average O/C ratio of clusters (column0) and sdev (column1) according to wave "+nameofwave(clusternumwave)
	Note/k HC_clst_avg "This wave stores the average H/C ratio of clusters (column0) and sdev (column1) according to wave "+nameofwave(clusternumwave)
	Note/k C_num_clst_avg "This wave stores the average Carbon number of clusters (column0) and sdev (column1) according to wave "+nameofwave(clusternumwave)
	Note/k N_num_clst_avg "This wave stores the average Nitrogen number of clusters (column0) and sdev (column1) according to wave "+nameofwave(clusternumwave)
	Note/k O_num_clst_avg "This wave stores the average Oxygen numberof clusters (column0) and sdev (column1) according to wave "+nameofwave(clusternumwave)
	Note/k H_num_clst_avg "This wave stores the average Hydrogen number of clusters (column0) and sdev (column1) according to wave "+nameofwave(clusternumwave)
	Note/k Mass_clst_avg "This wave stores the average absolute mass of clusters (column0) and sdev (column1) according to wave "+nameofwave(clusternumwave)
	Note/k Mass_clst_sum "This wave stores the sum of absolute mass of clusters according to wave "+nameofwave(clusternumwave)
	
End

Function AverageClstProperties_np(Clusternumwave, Clustersizewave)
	wave clusternumwave //-1 for unclustered ions and 1-x as cluster# for clustered ions
	wave clustersizewave //By default, the clustersize wave contains the number of unclustered ions in point 0

	//Assumes that the following waves already exists in the current folder
	wave Molarmass_cmp_fil
	wave N_num_fil
	wave C_num_fil
	wave H_num_fil
	wave O_num_fil
	wave OSc_fil
	wave OC_fil
	wave HC_fil
	wave Cmp_mass_fil
	
	if(waveexists(Molarmass_cmp_fil)==0 || waveexists(N_num_fil) ==0 || waveexists(C_num_fil) == 0|| waveexists(H_num_fil) == 0|| waveexists(O_num_fil) == 0|| waveexists(HC_fil) == 0|| waveexists(OC_fil) == 0|| waveexists(OSc_fil) == 0|| waveexists(Cmp_mass_fil) == 0)
		DoAlert 0, "Please make sure that waves Molarmass_cmp_fil, N_num_fil, C_num_fil, H_num_fil, O_num_fil, OSc_fil, OC_fil, HC_fil, Cmp_mass_fil are in the current folder"
		abort
	endif
	
	variable nrows = numpnts(clusternumwave)
	
	wavestats/q clusternumwave
	variable nclusters = V_max
	
	if(nclusters!=numpnts(clustersizewave)-1)
		print "Please make sure the clusternum wave and clustersize wave are results from the same clustering"
		abort
	endif
	
	variable i,j,k,clustercount
	variable numcount,tempsum
	
	make/o/d/n=(nclusters,2) Molarmass_clst_avg,OSc_clst_avg,OC_clst_avg,HC_clst_avg,C_num_clst_avg,N_num_clst_avg,O_num_clst_avg,H_num_clst_avg,Mass_clst_avg
	make/o/d/n=(nclusters) Mass_clst_sum
	
	for(i=0;i<nclusters;i+=1)
		make/o/d/n=(clustersizewave[i+1])/Free tempwave,tempwave1,tempwave2,tempwave3,tempwave4,tempwave5,tempwave6,tempwave7,tempwave8
		numcount=0
		for(j=0;j<nrows;j+=1)
			if(clusternumwave[j]==i+1)
				tempwave[numcount]=Molarmass_cmp_fil[j]
				tempwave1[numcount]=OSc_fil[j]
				tempwave2[numcount]=OC_fil[j]
				tempwave3[numcount]=HC_fil[j]
				tempwave4[numcount]=C_num_fil[j]
				tempwave5[numcount]=N_num_fil[j]
				tempwave6[numcount]=Cmp_mass_fil[j]
				tempwave7[numcount]=O_num_fil[j]
				tempwave8[numcount]=H_num_fil[j]
				numcount+=1
			endif
		endfor
		
			wavestats/q tempwave
			Molarmass_clst_avg[i][0]=V_avg
			Molarmass_clst_avg[i][1]=V_sdev
			wavestats/q tempwave1
			OSc_clst_avg[i][0]=V_avg
			OSc_clst_avg[i][1]=V_sdev
			wavestats/q tempwave2
			OC_clst_avg[i][0]=V_avg
			OC_clst_avg[i][1]=V_sdev
			wavestats/q tempwave3
			HC_clst_avg[i][0]=V_avg
			HC_clst_avg[i][1]=V_sdev
			wavestats/q tempwave4
			C_num_clst_avg[i][0]=V_avg
			C_num_clst_avg[i][1]=V_sdev
			wavestats/q tempwave5
			N_num_clst_avg[i][0]=V_avg
			N_num_clst_avg[i][1]=V_sdev
			wavestats/q tempwave6
			mass_clst_avg[i][0]=V_avg
			mass_clst_avg[i][1]=V_sdev
			Mass_clst_sum[i]=V_sum
			wavestats/q tempwave7
			O_num_clst_avg[i][0]=V_avg
			O_num_clst_avg[i][1]=V_sdev
			wavestats/q tempwave8
			H_num_clst_avg[i][0]=V_avg
			H_num_clst_avg[i][1]=V_sdev
	endfor
	
	Note/k Molarmass_clst_avg "This wave stores the average molarmass of clusters (column0) and sdev (column1) according to wave "+nameofwave(clusternumwave)
	Note/k OSc_clst_avg "This wave stores the average OSc of clusters (column0) and sdev (column1) according to wave "+nameofwave(clusternumwave)
	Note/k OC_clst_avg "This wave stores the average O/C ratio of clusters (column0) and sdev (column1) according to wave "+nameofwave(clusternumwave)
	Note/k HC_clst_avg "This wave stores the average H/C ratio of clusters (column0) and sdev (column1) according to wave "+nameofwave(clusternumwave)
	Note/k C_num_clst_avg "This wave stores the average Carbon number of clusters (column0) and sdev (column1) according to wave "+nameofwave(clusternumwave)
	Note/k N_num_clst_avg "This wave stores the average Nitrogen number of clusters (column0) and sdev (column1) according to wave "+nameofwave(clusternumwave)
	Note/k O_num_clst_avg "This wave stores the average Oxygen numberof clusters (column0) and sdev (column1) according to wave "+nameofwave(clusternumwave)
	Note/k H_num_clst_avg "This wave stores the average Hydrogen number of clusters (column0) and sdev (column1) according to wave "+nameofwave(clusternumwave)
	Note/k Mass_clst_avg "This wave stores the average absolute mass of clusters (column0) and sdev (column1) according to wave "+nameofwave(clusternumwave)
	Note/k Mass_clst_sum "This wave stores the sum of absolute mass of clusters according to wave "+nameofwave(clusternumwave)
	
End


//This function averages the selected properties for clusters (unweighted)
Function AverageClst_specific()

	string clusternumstr = "clusternumsort"
	string wave2avestr = "HC_fil"
	string newwavestr
	
	prompt clusternumstr, "Select wave of cluster#", popup WaveList("*",";","")
	prompt wave2avestr, "Select the property wave", popup WaveList("*",";","")
	prompt newwavestr, "Name of the new wave"
	Doprompt "Select waves",clusternumstr,wave2avestr,newwavestr
	
	if(V_flag)
		abort
	endif
	
	wave clusternum = $clusternumstr
	wave wave2ave = $wave2avestr
	
	wavestats clusternum
	variable nclsts=V_max
	variable npnts=numpnts(wave2ave)
	make/o/d/n=(nclsts) $newwavestr
	wave newwave=$newwavestr
	
	variable i,j,tempsum,count
	
	for(i=0;i<nclsts;i+=1)
		tempsum=0
		count=0
		for(j=0;j<npnts;j+=1)
			if(clusternum[j]==i+1)
				tempsum+=wave2ave[j]
				count+=1
			endif
		endfor
		newwave[i]=tempsum/count
	endfor

	Note/k newwave "This wave stores the unweighted averaged properties based on wave "+nameofwave(wave2ave)+" and cluster# according to wave "+nameofwave(clusternum)
	
End

Function AverageClst_specific_np(clusternum,wave2ave,newwavestr)
	wave clusternum,wave2ave
	string newwavestr
	
	wavestats clusternum
	variable nclsts=V_max
	variable npnts=numpnts(wave2ave)
	make/o/d/n=(nclsts) $newwavestr
	wave newwave=$newwavestr
	
	variable i,j,tempsum,count
	
	for(i=0;i<nclsts;i+=1)
		tempsum=0
		count=0
		for(j=0;j<npnts;j+=1)
			if(clusternum[j]==i+1)
				tempsum+=wave2ave[j]
				count+=1
			endif
		endfor
		newwave[i]=tempsum/count
	endfor

	Note/k newwave "This wave stores the unweighted averaged properties based on wave "+nameofwave(wave2ave)+" and cluster# according to wave "+nameofwave(clusternum)
	
End


//This function averages the selected properties for clusters (weighted according to weightwave)
Function AverageClst_specific_wt()
	
	string clusternumstr = "clusternumsort"
	string wave2avestr = "HC_fil"
	string weightwavestr = "cmp_mass_norm_fil"
	string newwavestr
	
	prompt clusternumstr, "Select wave of cluster#", popup WaveList("*",";","")
	prompt wave2avestr, "Select the property wave", popup WaveList("*",";","")
	prompt weightwavestr, "Select the weighting wave", popup WaveList("*",";","")
	prompt newwavestr, "Name of the new wave"
	Doprompt "Select waves",clusternumstr,wave2avestr,weightwavestr, newwavestr
	
	if(V_flag)
		abort
	endif
	
	wave clusternum = $clusternumstr
	wave wave2ave = $wave2avestr
	wave weightwave = $weightwavestr
	
	wavestats clusternum
	variable nclsts=V_max
	variable npnts=numpnts(wave2ave)
	make/o/d/n=(nclsts) $newwavestr
	wave newwave=$newwavestr
	
	variable i,j,tempsum,tempweightsum,count
	
	for(i=0;i<nclsts;i+=1)
		tempsum=0
		tempweightsum=0
		count=0
		for(j=0;j<npnts;j+=1)
			if(clusternum[j]==i+1)
				tempsum+=wave2ave[j]*weightwave[j]
				tempweightsum+=weightwave[j]
				count+=1
			endif
		endfor
		newwave[i]=tempsum/tempweightsum
	endfor

	Note/k newwave "This wave stores the weighted averaged properties based on wave "+nameofwave(wave2ave)+" , weighting based on wave "+nameofwave(weightwave)+" and cluster# according to wave "+nameofwave(clusternum)
	
End

Function AverageClst_specific_wt_np(clusternum,wave2ave,weightwave,newwavestr)
	wave clusternum,wave2ave,weightwave
	string newwavestr
	
	wavestats clusternum
	variable nclsts=V_max
	variable npnts=numpnts(wave2ave)
	make/o/d/n=(nclsts) $newwavestr
	wave newwave=$newwavestr
	
	variable i,j,tempsum,tempweightsum,count
	
	for(i=0;i<nclsts;i+=1)
		tempsum=0
		tempweightsum=0
		count=0
		for(j=0;j<npnts;j+=1)
			if(clusternum[j]==i+1)
				tempsum+=wave2ave[j]*weightwave[j]
				tempweightsum+=weightwave[j]
				count+=1
			endif
		endfor
		newwave[i]=tempsum/tempweightsum
	endfor

	Note/k newwave "This wave stores the weighted averaged properties based on wave "+nameofwave(wave2ave)+" , weighting based on wave "+nameofwave(weightwave)+" and cluster# according to wave "+nameofwave(clusternum)
	
End

//This function separates mass spectrumfor each cluster and stores spectra of clusters in columns of a matrix
Function ClusterSpectraMatrix()
	string signalwavestr = "cmp_mass_norm_fil"
	string clusternumstr = "clusternumsort"
	string suffix = ""
	
	prompt signalwavestr, "Select wave of signal/mass", popup WaveList("*",";","")
	prompt clusternumstr, "Select wave of cluster#", popup WaveList("*",";","")
	prompt suffix, "Suffix to add to the new matrix"
	Doprompt "Select waves",signalwavestr,clusternumstr,suffix
	
	if(V_flag)
		abort
	endif

	wave signalwave = $signalwavestr
	wave clusternum = $clusternumstr
	
	variable ncols=numpnts(signalwave)
	wavestats/q clusternum
	variable nclusters=V_max
	
	if(ncols!=numpnts(clusternum))
		DoAlert 0, "Please make sure wave clusternum and signalwave have the same length"
		abort
	endif
	
	make/o/d/n=(nclusters,ncols) $("ClusterSpectra"+suffix)
	wave clusterspectra=$("ClusterSpectra"+suffix)
	clusterspectra=nan
	
	variable i,j
	
	for(i=0;i<ncols;i+=1)
		if(clusternum[i]!=-1)
			clusterspectra[clusternum[i]-1][i]=signalwave[i]
		endif
	endfor
	
	Note/k clusterspectra "This matrix stores the spectrum of each cluster in each column"
	Note clusterspectra "Function used: ClusterSpectraMatrix. signalwave <- "+nameofwave(signalwave)+" ; clusternum <- "+nameofwave(clusternum)
	
End


Function ClusterSpectraMatrix_np(signalwave,clusternum,suffix)
	wave signalwave
	wave clusternum
	string suffix
	
	variable ncols=numpnts(signalwave)
	wavestats/q clusternum
	variable nclusters=V_max
	
	if(ncols!=numpnts(clusternum))
		DoAlert 0, "Please make sure wave clusternum and signalwave have the same length"
		abort
	endif
	
	make/o/d/n=(nclusters,ncols) $("ClusterSpectra"+suffix)
	wave clusterspectra=$("ClusterSpectra"+suffix)
	clusterspectra=nan
	
	variable i,j
	
	for(i=0;i<ncols;i+=1)
		if(clusternum[i]!=-1)
			clusterspectra[clusternum[i]-1][i]=signalwave[i]
		endif
	endfor
	
	Note/k clusterspectra "This matrix stores the spectrum of each cluster in each column"
	Note clusterspectra "Function used: ClusterSpectraMatrix_np. signalwave <- "+nameofwave(signalwave)+" ; clusternum <- "+nameofwave(clusternum)
	
End

//This function graphs the four parameters that helps determine the optimal eps 
Function GraphOptEpsParas()
	
	string Nc_allstr = "NumofClusters_all"
	string Nc_no1Mstr = "NumofClusters_no1M"
	string UngroupedMstr = "ungroupedMass"
	string R_interSDstr = "Ratio_avginterSD_eps"
	string epsstr = "epsilon"
	
	prompt Nc_allstr, "Select wave of Nc (all)", popup WaveList("*",";","")
	prompt Nc_no1Mstr, "Select wave of Nc (w/o one-member)", popup WaveList("*",";","")
	prompt ungroupedMstr, "Select wave of unclustered mass", popup WaveList("*",";","")
	prompt R_interSDstr, "Select wave of R_interED", popup WaveList("*",";","")
	prompt epsstr, "Select wave of epsilon", popup WaveList("*",";","")
	Doprompt "Select waves",Nc_allstr,Nc_no1Mstr,UngroupedMstr,R_interSDstr,epsstr
	
	if(V_flag)
		abort
	endif
	
	wave Nc_all = $Nc_allstr
	wave Nc_no1M = $Nc_no1Mstr
	wave UngroupedM = $ungroupedMstr
	wave R_interSD = $R_interSDstr
	wave eps = $epsstr

	display/l=Nc Nc_all vs eps
	appendtograph/l=Nc Nc_no1M vs eps
	appendtograph/l=fm UngroupedM vs eps
	appendtograph/r=Rinter R_interSD vs eps
	
	SetAxis Nc 0,35
	SetAxis fm 0,0.2
	SetAxis Rinter 0,4
	
	Label Nc "\\f02N\\Bc"
	Label fm "\\f02f\\Bm,unclustered"
	Label Rinter "\\f02R\\BinterClst"
	Label bottom "\\Z16\\F'Symbol'e"
	
	ModifyGraph mode=4,msize=4,lsize=1.5,marker($Nc_allstr)=19
	ModifyGraph marker($Nc_no1Mstr)=8,marker($ungroupedMstr)=17
	ModifyGraph rgb($ungroupedMstr)=(65535,43690,0)
	ModifyGraph marker($R_interSDstr)=29
	ModifyGraph rgb($R_interSDstr)=(0,0,65535)

	ModifyGraph fSize=12,axRGB(Nc)=(65535,0,0),axRGB(fm)=(65535,43690,0)
	ModifyGraph axRGB(Rinter)=(0,0,65535),tlblRGB(Nc)=(65535,0,0)
	ModifyGraph tlblRGB(fm)=(65535,43690,0),tlblRGB(Rinter)=(0,0,65535)
	ModifyGraph alblRGB(Nc)=(65535,0,0),alblRGB(fm)=(65535,43690,0)
	ModifyGraph alblRGB(Rinter)=(0,0,65535)
	
	ModifyGraph lblPosMode(Nc)=4,lblPosMode(fm)=4,lblPosMode(Rinter)=4
	ModifyGraph lblPos(Nc)=50,lblPos(fm)=50,lblPos(Rinter)=50
	ModifyGraph axisEnab(bottom)={0.2,0.95}
	ModifyGraph freePos(Nc)={0.05,kwFraction}
	ModifyGraph freePos(fm)={0.2,kwFraction}
	ModifyGraph freePos(Rinter)={0.05,kwFraction}
	
	Legend/C/N=text0/J "\\s("+Nc_Allstr+") N\\Bc \\Mwith one-member\r\\s("+Nc_no1Mstr+") N\\Bc \\Mwithout one-member"
	AppendText "\\s("+ungroupedMstr+") Unclustered mass fraction\r\\s("+R_interSDstr+") AvginterED/\\F'Symbol'e"

	modifygraph width=450, height=250

End

Function GraphOptEpsParas_np(Nc_all,Nc_no1M,ungroupedM,R_interSD,eps)
	wave Nc_all
	wave Nc_no1M
	wave UngroupedM
	wave R_interSD
	wave eps
	
	string Nc_allstr = nameofwave(Nc_all)
	string Nc_no1Mstr = nameofwave(Nc_no1M)
	string UngroupedMstr = nameofwave(ungroupedM)
	string R_interSDstr = nameofwave(R_interSD)
	

	display/l=Nc Nc_all vs eps
	appendtograph/l=Nc Nc_no1M vs eps
	appendtograph/l=fm UngroupedM vs eps
	appendtograph/r=Rinter R_interSD vs eps
	
	SetAxis Nc 0,35
	SetAxis fm 0,0.2
	SetAxis Rinter 0,4
	
	Label Nc "\\f02N\\Bc"
	Label fm "\\f02f\\Bm,unclustered"
	Label Rinter "\\f02R\\BinterClst"
	Label bottom "\\Z16\\F'Symbol'e"
	
	ModifyGraph mode=4,msize=4,lsize=1.5,marker($Nc_allstr)=19
	ModifyGraph marker($Nc_no1Mstr)=8,marker($ungroupedMstr)=17
	ModifyGraph rgb($ungroupedMstr)=(65535,43690,0)
	ModifyGraph marker($R_interSDstr)=29
	ModifyGraph rgb($R_interSDstr)=(0,0,65535)

	ModifyGraph fSize=12,axRGB(Nc)=(65535,0,0),axRGB(fm)=(65535,43690,0)
	ModifyGraph axRGB(Rinter)=(0,0,65535),tlblRGB(Nc)=(65535,0,0)
	ModifyGraph tlblRGB(fm)=(65535,43690,0),tlblRGB(Rinter)=(0,0,65535)
	ModifyGraph alblRGB(Nc)=(65535,0,0),alblRGB(fm)=(65535,43690,0)
	ModifyGraph alblRGB(Rinter)=(0,0,65535)
	
	ModifyGraph lblPosMode(Nc)=4,lblPosMode(fm)=4,lblPosMode(Rinter)=4
	ModifyGraph lblPos(Nc)=50,lblPos(fm)=50,lblPos(Rinter)=50
	ModifyGraph axisEnab(bottom)={0.2,0.95}
	ModifyGraph freePos(Nc)={0.05,kwFraction}
	ModifyGraph freePos(fm)={0.2,kwFraction}
	ModifyGraph freePos(Rinter)={0.05,kwFraction}
	
	Legend/C/N=text0/J "\\s("+Nc_Allstr+") N\\Bc \\Mwith one-member\r\\s("+Nc_no1Mstr+") N\\Bc \\Mwithout one-member"
	AppendText "\\s("+ungroupedMstr+") Unclustered mass fraction\r\\s("+R_interSDstr+") AvginterED/\\F'Symbol'e"

	modifygraph width=450, height=250

End


//This function calculate the thermograms of clusters by adding up th thermograms of absolute signals of each cluster
//The unclustered ions are not considered
Function SumThermogramClst()
	
	string ymatrixstr = "cmp_matrix_fil"
	string clusternumstr = "clusternumsort"
	string suffix
	
	prompt ymatrixstr, "Select matrix of absolute signals", popup WaveList("*",";","")
	prompt clusternumstr, "Select wave of cluster#", popup WaveList("*",";","")
	prompt suffix, "Suffix to add to the new matrix"
	Doprompt "Select waves",ymatrixstr, clusternumstr, suffix
	
	if(V_flag)
		abort
	endif
	
	wave ymatrix =$ymatrixstr 
	wave clusternum = $clusternumstr
	
	variable ncols=dimsize(ymatrix,1)
	variable nrows=dimsize(ymatrix,0)
	if(ncols!=numpnts(clusternum))
		print "Please make sure the number of columns of the matrix equals to the number of rows of the wave "
		abort
	endif
	
	wavestats/q clusternum
	make/o/d/n=(nrows,V_max) $("SumTG_clst"+suffix)
	wave sumthermograms=$("SumTG_clst"+suffix)
	sumthermograms=0
	
	variable i,j
	
	for(j=0;j<ncols;j+=1)
		
		if(clusternum[j]!=-1)	
			sumthermograms[][clusternum[j]-1]+=ymatrix[p][j]
		endif
	
	endfor
	
	Note/k sumthermograms "Each column in this matrix stores the sum of TG with abs signals/mass of a cluster"
	Note sumthermograms "Function used: SumThermogramClst. ymatrix <- "+ymatrixstr+" and Clusternum <- "+clusternumstr
	
end

Function SumThermogramClst_np(ymatrix,clusternum,suffix)
	wave ymatrix //matrix with absolute signals/mass
	wave clusternum
	string suffix
	
	variable ncols=dimsize(ymatrix,1)
	variable nrows=dimsize(ymatrix,0)
	if(ncols!=numpnts(clusternum))
		print "Please make sure the number of columns of the matrix equals to the number of rows of the wave "
		abort
	endif
	
	wavestats/q clusternum
	make/o/d/n=(nrows,V_max) $("SumTG_clst"+suffix)
	wave sumthermograms=$("SumTG_clst"+suffix)
	sumthermograms=0
	
	variable i,j
	
	for(j=0;j<ncols;j+=1)
		
		if(clusternum[j]!=-1)	
			sumthermograms[][clusternum[j]-1]+=ymatrix[p][j]
		endif
	
	endfor
	
	Note/k sumthermograms "Each column in this matrix stores the sum of TG with abs signals/mass of a cluster"
	Note sumthermograms "Function used: SumThermogramClst_np. ymatrix <- "+nameofwave(ymatrix)+" and Clusternum <- "+nameofwave(clusternum)
	
	
end

//This function counts cluster# starting with 0, instead of 1
Function SumThermogramClst0_np(ymatrix,clusternum,suffix)
	wave ymatrix //matrix with absolute signals/mass
	wave clusternum
	string suffix
	
	variable ncols=dimsize(ymatrix,1)
	variable nrows=dimsize(ymatrix,0)
	if(ncols!=numpnts(clusternum))
		print "Please make sure the number of columns of the matrix equals to the number of rows of the wave "
		abort
	endif
	
	wavestats/q clusternum
	make/o/d/n=(nrows,V_max+1) $("SumTG_clst"+suffix)
	wave sumthermograms=$("SumTG_clst"+suffix)
	sumthermograms=0
	
	variable i,j
	
	for(j=0;j<ncols;j+=1)
		
		if(clusternum[j]>=0)	
			sumthermograms[][clusternum[j]]+=ymatrix[p][j]
		endif
	
	endfor
	
	Note/k sumthermograms "Each column in this matrix stores the sum of TG with abs signals/mass of a cluster"
	Note sumthermograms "Function used: SumThermogramClst0_np. ymatrix <- "+nameofwave(ymatrix)+" and Clusternum <- "+nameofwave(clusternum)
	
	
end



//This function plots cluster 1-x, plots both weighted and unweighted average thermograms on the same figure
//Each cluster ends up in a different graph
Function GraphClstTG2_MultiGraphs()
	
	string clusternumstr = "Clusternumsort"
	string ymatrixstr = "cmp_matrix_sn_10a_fil"
	string ymatrix_avestr = "Clusteravessort"
	string ymatrix_ave_wtstr = "Clusterave_wtsort"
	string xwavestr = "temp_c_10a"
	
	prompt clusternumstr, "Select wave of cluster#", popup WaveList("*",";","")
	prompt ymatrixstr, "Select matrix of individual TGs", popup WaveList("*",";","")
	prompt ymatrix_avestr, "Select matrix of avg TG of clusters", popup WaveList("*",";","")
	prompt ymatrix_ave_wtstr, "Select matrix of weighted avg TG", popup WaveList("*",";","")
	prompt xwavestr, "Select temperature wave", popup WaveList("*",";","")
	Doprompt "Select waves",clusternumstr, ymatrixstr, ymatrix_avestr, ymatrix_ave_wtstr, xwavestr
	
	if(V_flag)
		abort
	endif
	
	wave clusternum = $clusternumstr
	wave ymatrix = $ymatrixstr
	wave ymatrix_ave = $ymatrix_avestr
	wave ymatrix_ave_wt = $ymatrix_ave_wtstr
	wave xwave = $xwavestr
	
	variable nrows=dimsize(ymatrix,0)
	variable ncols=dimsize(ymatrix,1)
	wavestats/q clusternum
	variable nclusters=V_max
	 variable i,j,k
	
	colortab2wave dBZ21
	wave M_colors
	SetScale/I x,1,nclusters,M_colors
	variable red,green,blue
	variable indexpnt
	variable zlevel
	
	for(i=0;i<nclusters;i+=1)
		k=0
		zlevel=i+1
		for(j=0;j<ncols;j+=1)
			if(ClusterNum[j]==i+1)
				if(k==0)
					string graphname = "Cluster "+num2istr(i+1)
					display/N=$graphname as graphname
					appendtograph ymatrix[][j] vs xwave
				else
					appendtograph ymatrix[][j] vs xwave
				endif
				k+=1
			endif
		endfor
		red=M_colors(zlevel)[0]
		green=M_colors(zlevel)[1]
		blue=M_colors(zlevel)[2]
		ModifyGraph rgb=(red,green,blue)
		appendtograph ymatrix_ave[][i] vs xwave
		appendtograph ymatrix_Ave_wt[][i] vs xwave
		ModifyGraph lsize($(nameofwave(ymatrix_ave)))=2.5, lsize($(nameofwave(ymatrix_ave_wt)))=2.5, rgb($(nameofwave(ymatrix_ave_wt)))=(0,0,0),rgb($(nameofwave(ymatrix_ave)))=(30583,30583,30583)
		SetAxis left -0.1,1.3
		SetAxis bottom 30,205
		ModifyGraph mirror=1
	endfor
End


Function GraphClstTG2_MultiGraphs_np(Clusternum,ymatrix,ymatrix_ave,ymatrix_ave_wt,xwave)
	wave clusternum
	wave ymatrix
	wave ymatrix_ave
	wave ymatrix_ave_wt
	wave xwave
	
	variable nrows=dimsize(ymatrix,0)
	variable ncols=dimsize(ymatrix,1)
	wavestats/q clusternum
	variable nclusters=V_max
	 variable i,j,k
	
	colortab2wave dBZ21
	wave M_colors
	SetScale/I x,1,nclusters,M_colors
	variable red,green,blue
	variable indexpnt
	variable zlevel
	
	for(i=0;i<nclusters;i+=1)
		k=0
		zlevel=i+1
		for(j=0;j<ncols;j+=1)
			if(ClusterNum[j]==i+1)
				if(k==0)
					string graphname = "Cluster "+num2istr(i+1)
					display/N=$graphname as graphname
					appendtograph ymatrix[][j] vs xwave
				else
					appendtograph ymatrix[][j] vs xwave
				endif
				k+=1
			endif
		endfor
		red=M_colors(zlevel)[0]
		green=M_colors(zlevel)[1]
		blue=M_colors(zlevel)[2]
		ModifyGraph rgb=(red,green,blue)
		appendtograph ymatrix_ave[][i] vs xwave
		appendtograph ymatrix_Ave_wt[][i] vs xwave
		ModifyGraph lsize($(nameofwave(ymatrix_ave)))=2.5, lsize($(nameofwave(ymatrix_ave_wt)))=2.5, rgb($(nameofwave(ymatrix_ave_wt)))=(0,0,0),rgb($(nameofwave(ymatrix_ave)))=(30583,30583,30583)
		SetAxis left -0.1,1.3
		SetAxis bottom 30,205
		ModifyGraph mirror=1
	endfor
End


//This function plots cluster 1-x, plots both weighted and unweighted average thermograms, and the individual thermograms on the same figure
//Each cluster ends up as a different panel in the same graph 
Function GraphClstTG2_OneGraph()

	string clusternumstr = "Clusternumsort"
	string ymatrixstr = "cmp_matrix_sn_10a_fil"
	string ymatrix_avestr = "Clusterave_unwtsort"
	string ymatrix_ave_wtstr = "Clusterave_wtsort"
	string xwavestr = "temp_c_10a"
	variable ncols = 3
	
	prompt clusternumstr, "Select wave of cluster#", popup WaveList("*",";","")
	prompt ymatrixstr, "Select matrix of individual TGs", popup WaveList("*",";","")
	prompt ymatrix_avestr, "Select matrix of avg TG of clusters", popup WaveList("*",";","")
	prompt ymatrix_ave_wtstr, "Select matrix of weighted avg TG", popup WaveList("*",";","")
	prompt xwavestr, "Select temperature wave", popup WaveList("*",";","")
	prompt ncols, "Columns of panels of graph"
	Doprompt "Select waves",clusternumstr, ymatrixstr, ymatrix_avestr, ymatrix_ave_wtstr, xwavestr, ncols
	
	if(V_flag)
		abort
	endif
	
	wave clusternum = $clusternumstr
	wave ymatrix = $ymatrixstr
	wave ymatrix_ave = $ymatrix_avestr
	wave ymatrix_ave_wt = $ymatrix_ave_wtstr
	wave xwave = $xwavestr

	variable nclsts = dimsize(ymatrix_ave,1) // number of clusters
	variable nions = dimsize(ymatrix,1) //number of individual ions
	variable nrows = ceil(nclsts/ncols) // number of rows
	
	colortab2wave dBZ21
	wave M_colors
	SetScale/I x,1,nclsts,M_colors
	variable red,green,blue
	
	variable col_offset
	variable col_len 
	variable col_space 
	variable col_multiplier = 2
	
	col_offset = 0.05
	//Let col_len=col_multiplier*col_space, then col_offset+ncols*col_len+(ncols-1)*col_space=1
	//so that
	col_space = (1-col_offset)/((1+col_multiplier)*ncols-1)
	col_len=col_multiplier*col_space
	
	variable row_len
	variable row_space = 0.2/nrows
	row_len = 1/nrows
	
	string yaxstr, xaxstr
	variable i, j,k
	variable wvcount = 0
	variable ioncount = 0
	variable index
	
	variable xlowFra,xhighFra,ylowFra,yhighFra
	string abcstr="a;b;c;d;e;f;g;h;i;j;k;l;m;n;o;p;q;r;s;t;u;v;w;x;y;z"
	string NameofPanel
	
	string graphname = "AvgThermogramsofClst"
	
	display/N=$graphname as graphname
	modifygraph width=500, height=500

	for(i=0;i<nrows;i+=1)
		for(j=0;j<ncols;j+=1)
			if(wvcount > nclsts-1)
				return 1
			endif
			index=j+1+ncols*i
			yaxstr = "y" + num2istr(index)
			xaxstr = "x" + num2istr(index)
			red=M_colors(index)[0]
			green=M_colors(index)[1]
			blue=M_colors(index)[2]
			for(k=0;k<nions;k+=1)
				if(clusternum[k]==index)
					appendtograph/l=$yaxstr/b=$xaxstr ymatrix[][k] vs xwave
					if(ioncount == 0)
						ModifyGraph rgb($(ymatrixstr))=(red,green,blue)
					else
						ModifyGraph rgb($(ymatrixstr+"#"+num2istr(ioncount)))=(red,green,blue)
					endif
					ioncount+=1
				endif
			endfor
			appendtograph/l=$yaxstr/b=$xaxstr ymatrix_ave[][index-1] vs xwave
			appendtograph/l=$yaxstr/b=$xaxstr ymatrix_Ave_wt[][index-1] vs xwave
			if(index==1)
				ModifyGraph lsize($(ymatrix_avestr))=2.5, lsize($(ymatrix_ave_wtstr))=2.5, rgb($(ymatrix_ave_wtstr))=(0,0,0),rgb($(ymatrix_avestr))=(30583,30583,30583)
			else
				ModifyGraph lsize($(ymatrix_avestr+"#"+num2istr(index-1)))=2.5, lsize($(ymatrix_ave_wtstr+"#"+num2istr(index-1)))=2.5, rgb($(ymatrix_ave_wtstr+"#"+num2istr(index-1)))=(0,0,0),rgb($(ymatrix_avestr+"#"+num2istr(index-1)))=(30583,30583,30583)
			endif
			if(j==0)
				//ModifyGraph axisEnab($yaxstr)={1-row_len*(i+1)+row_space*(nrows-i),1-(row_len)*i-row_space*(i)},axisEnab($xaxstr)={col_offset+col_len*j,col_offset+col_len*(j+1)} // col_len*j+col_space*j,col_len2*(j+1)-col_space*(ncols-j-1)}
				ylowFra = 1-row_len*(i+1)+row_space //*(nrows-i)
				yhighFra = 1-(row_len)*i-row_space //*(i)
				xlowFra = col_offset+col_len*j
				xhighFra = col_offset+col_len*(j+1)
				ModifyGraph axisEnab($yaxstr)={ylowFra,yhighFra},axisEnab($xaxstr)={xlowFra,xhighFra}
			else
				ylowFra = 1-row_len*(i+1)+row_space //*(nrows-i)
				yhighFra = 1-(row_len)*i-row_space //*(i)
				xlowFra = col_offset+col_len*j+col_space*j
				xhighFra = col_offset+col_len*(j+1)+col_space*j
				ModifyGraph axisEnab($yaxstr)={ylowFra,yhighFra},axisEnab($xaxstr)={xlowFra,xhighFra}
			endif
			//print ylowFra,yhighFra,xlowFra,xhighFra
			SetAxis $yaxstr -0.1,1.3
			SetAxis $xaxstr 30,205			
			ModifyGraph freePos($yaxstr)={0,$xaxstr},freePos($xaxstr)={0,$yaxstr}
			//ModifyGraph noLabel($xaxstr)=1//, noLabel($yaxstr)=1
			//NameofPanel= "\\Z12"+stringfromlist(index-1,abcstr,";")+")"
			NameofPanel= "\\Z12Clst#"+num2istr(index)
			TextBox/C/F=0/B=1/N=$("name"+num2istr(index))/A=MC/X=((xhighFra-col_space/2)*100-50)/Y=(yhighFra*100-50) nameofpanel
			wvcount +=1
		endfor
	endfor

	
End

Function GraphClstTG2_OneGraph_np(clusternum, ymatrix,ymatrix_ave,ymatrix_ave_wt, xwave, ncols)

	wave clusternum
	wave ymatrix
	wave ymatrix_ave 
	wave ymatrix_ave_wt 
	wave xwave
	variable ncols

	string clusternumstr = nameofwave(clusternum)
	string ymatrixstr = nameofwave(ymatrix)
	string ymatrix_avestr = nameofwave(ymatrix_ave)
	string ymatrix_ave_wtstr = nameofwave(ymatrix_ave_wt)
	string xwavestr = nameofwave(xwave)

	variable nclsts = dimsize(ymatrix_ave,1) // number of clusters
	variable nions = dimsize(ymatrix,1) //number of individual ions
	variable nrows = ceil(nclsts/ncols) // number of rows
	
	colortab2wave dBZ21
	wave M_colors
	SetScale/I x,1,nclsts,M_colors
	variable red,green,blue
	
	variable col_offset
	variable col_len 
	variable col_space 
	variable col_multiplier = 2
	
	col_offset = 0.05
	//Let col_len=col_multiplier*col_space, then col_offset+ncols*col_len+(ncols-1)*col_space=1
	//so that
	col_space = (1-col_offset)/((1+col_multiplier)*ncols-1)
	col_len=col_multiplier*col_space
	
	variable row_len
	variable row_space = 0.2/nrows
	row_len = 1/nrows
	
	string yaxstr, xaxstr
	variable i, j,k
	variable wvcount = 0
	variable ioncount = 0
	variable index
	
	variable xlowFra,xhighFra,ylowFra,yhighFra
	string abcstr="a;b;c;d;e;f;g;h;i;j;k;l;m;n;o;p;q;r;s;t;u;v;w;x;y;z"
	string NameofPanel
	
	string graphname = "AvgThermogramsofClst"
	
	display/N=$graphname as graphname
	modifygraph width=500, height=500

	for(i=0;i<nrows;i+=1)
		for(j=0;j<ncols;j+=1)
			if(wvcount > nclsts-1)
				return 1
			endif
			index=j+1+ncols*i
			yaxstr = "y" + num2istr(index)
			xaxstr = "x" + num2istr(index)
			red=M_colors(index)[0]
			green=M_colors(index)[1]
			blue=M_colors(index)[2]
			for(k=0;k<nions;k+=1)
				if(clusternum[k]==index)
					appendtograph/l=$yaxstr/b=$xaxstr ymatrix[][k] vs xwave
					if(ioncount == 0)
						ModifyGraph rgb($(ymatrixstr))=(red,green,blue)
					else
						ModifyGraph rgb($(ymatrixstr+"#"+num2istr(ioncount)))=(red,green,blue)
					endif
					ioncount+=1
				endif
			endfor
			appendtograph/l=$yaxstr/b=$xaxstr ymatrix_ave[][index-1] vs xwave
			appendtograph/l=$yaxstr/b=$xaxstr ymatrix_Ave_wt[][index-1] vs xwave
			if(index==1)
				ModifyGraph lsize($(ymatrix_avestr))=2.5, lsize($(ymatrix_ave_wtstr))=2.5, rgb($(ymatrix_ave_wtstr))=(0,0,0),rgb($(ymatrix_avestr))=(30583,30583,30583)
			else
				ModifyGraph lsize($(ymatrix_avestr+"#"+num2istr(index-1)))=2.5, lsize($(ymatrix_ave_wtstr+"#"+num2istr(index-1)))=2.5, rgb($(ymatrix_ave_wtstr+"#"+num2istr(index-1)))=(0,0,0),rgb($(ymatrix_avestr+"#"+num2istr(index-1)))=(30583,30583,30583)
			endif
			if(j==0)
				//ModifyGraph axisEnab($yaxstr)={1-row_len*(i+1)+row_space*(nrows-i),1-(row_len)*i-row_space*(i)},axisEnab($xaxstr)={col_offset+col_len*j,col_offset+col_len*(j+1)} // col_len*j+col_space*j,col_len2*(j+1)-col_space*(ncols-j-1)}
				ylowFra = 1-row_len*(i+1)+row_space //*(nrows-i)
				yhighFra = 1-(row_len)*i-row_space //*(i)
				xlowFra = col_offset+col_len*j
				xhighFra = col_offset+col_len*(j+1)
				ModifyGraph axisEnab($yaxstr)={ylowFra,yhighFra},axisEnab($xaxstr)={xlowFra,xhighFra}
			else
				ylowFra = 1-row_len*(i+1)+row_space //*(nrows-i)
				yhighFra = 1-(row_len)*i-row_space //*(i)
				xlowFra = col_offset+col_len*j+col_space*j
				xhighFra = col_offset+col_len*(j+1)+col_space*j
				ModifyGraph axisEnab($yaxstr)={ylowFra,yhighFra},axisEnab($xaxstr)={xlowFra,xhighFra}
			endif
			//print ylowFra,yhighFra,xlowFra,xhighFra
			SetAxis $yaxstr -0.1,1.3
			SetAxis $xaxstr 30,205			
			ModifyGraph freePos($yaxstr)={0,$xaxstr},freePos($xaxstr)={0,$yaxstr}
			//ModifyGraph noLabel($xaxstr)=1//, noLabel($yaxstr)=1
			//NameofPanel= "\\Z12"+stringfromlist(index-1,abcstr,";")+")"
			NameofPanel= "\\Z12Clst#"+num2istr(index)
			TextBox/C/F=0/B=1/N=$("name"+num2istr(index))/A=MC/X=((xhighFra-col_space/2)*100-50)/Y=(yhighFra*100-50) nameofpanel
			wvcount +=1
		endfor
	endfor

	
End

//This function plots cluster 1-x, plots the average and the individual thermograms on the same figure
//Each cluster ends up as a different panel in the same graph 
Function GraphClstTG1_OneGraph()

	string clusternumstr = "Clusternumsort"
	string ymatrixstr = "cmp_matrix_sn_10a_fil"
	string ymatrix_avestr = "Clusterave_unwt"
	string xwavestr = "temp_c_10a"
	variable ncols = 3
	
	prompt clusternumstr, "Select wave of cluster#", popup WaveList("*",";","")
	prompt ymatrixstr, "Select matrix of individual TGs", popup WaveList("*",";","")
	prompt ymatrix_avestr, "Select matrix of avg TG of clusters", popup WaveList("*",";","")
	prompt xwavestr, "Select temperature wave", popup WaveList("*",";","")
	prompt ncols, "Columns of panels of graph"
	Doprompt "Select waves",clusternumstr, ymatrixstr, ymatrix_avestr, xwavestr, ncols
	
	if(V_flag)
		abort
	endif
	
	wave clusternum = $clusternumstr
	wave ymatrix = $ymatrixstr
	wave ymatrix_ave = $ymatrix_avestr
	wave xwave = $xwavestr

	variable nclsts = dimsize(ymatrix_ave,1) // number of clusters
	variable nions = dimsize(ymatrix,1) //number of individual ions
	variable nrows = ceil(nclsts/ncols) // number of rows
	
	colortab2wave dBZ21
	wave M_colors
	SetScale/I x,1,nclsts,M_colors
	variable red,green,blue
	
	variable col_offset
	variable col_len 
	variable col_space 
	variable col_multiplier = 2
	
	col_offset = 0.05
	//Let col_len=col_multiplier*col_space, then col_offset+ncols*col_len+(ncols-1)*col_space=1
	//so that
	col_space = (1-col_offset)/((1+col_multiplier)*ncols-1)
	col_len=col_multiplier*col_space
	
	variable row_len
	variable row_space = 0.2/nrows
	row_len = 1/nrows
	
	string yaxstr, xaxstr
	variable i, j,k
	variable wvcount = 0
	variable ioncount = 0
	variable index
	
	variable xlowFra,xhighFra,ylowFra,yhighFra
	string abcstr="a;b;c;d;e;f;g;h;i;j;k;l;m;n;o;p;q;r;s;t;u;v;w;x;y;z"
	string NameofPanel
	
	string graphname = "AvgThermogramsofClst"
	
	display/N=$graphname as graphname
	modifygraph width=500, height=500

	for(i=0;i<nrows;i+=1)
		for(j=0;j<ncols;j+=1)
			if(wvcount > nclsts-1)
				return 1
			endif
			index=j+1+ncols*i
			yaxstr = "y" + num2istr(index)
			xaxstr = "x" + num2istr(index)
			red=M_colors(index)[0]
			green=M_colors(index)[1]
			blue=M_colors(index)[2]
			for(k=0;k<nions;k+=1)
				if(clusternum[k]==index)
					appendtograph/l=$yaxstr/b=$xaxstr ymatrix[][k] vs xwave
					if(ioncount == 0)
						ModifyGraph rgb($(ymatrixstr))=(red,green,blue)
					else
						ModifyGraph rgb($(ymatrixstr+"#"+num2istr(ioncount)))=(red,green,blue)
					endif
					ioncount+=1
				endif
			endfor
			appendtograph/l=$yaxstr/b=$xaxstr ymatrix_ave[][index-1] vs xwave
			if(index==1)
				ModifyGraph lsize($(ymatrix_avestr))=2.5, rgb($(ymatrix_avestr))=(0,0,0)
			else
				ModifyGraph lsize($(ymatrix_avestr+"#"+num2istr(index-1)))=2.5, rgb($(ymatrix_avestr+"#"+num2istr(index-1)))=(0,0,0)
			endif
			if(j==0)
				//ModifyGraph axisEnab($yaxstr)={1-row_len*(i+1)+row_space*(nrows-i),1-(row_len)*i-row_space*(i)},axisEnab($xaxstr)={col_offset+col_len*j,col_offset+col_len*(j+1)} // col_len*j+col_space*j,col_len2*(j+1)-col_space*(ncols-j-1)}
				ylowFra = 1-row_len*(i+1)+row_space //*(nrows-i)
				yhighFra = 1-(row_len)*i-row_space //*(i)
				xlowFra = col_offset+col_len*j
				xhighFra = col_offset+col_len*(j+1)
				ModifyGraph axisEnab($yaxstr)={ylowFra,yhighFra},axisEnab($xaxstr)={xlowFra,xhighFra}
			else
				ylowFra = 1-row_len*(i+1)+row_space //*(nrows-i)
				yhighFra = 1-(row_len)*i-row_space //*(i)
				xlowFra = col_offset+col_len*j+col_space*j
				xhighFra = col_offset+col_len*(j+1)+col_space*j
				ModifyGraph axisEnab($yaxstr)={ylowFra,yhighFra},axisEnab($xaxstr)={xlowFra,xhighFra}
			endif
			//print ylowFra,yhighFra,xlowFra,xhighFra
			SetAxis $yaxstr -0.1,1.3
			SetAxis $xaxstr 30,205			
			ModifyGraph freePos($yaxstr)={0,$xaxstr},freePos($xaxstr)={0,$yaxstr}
			//ModifyGraph noLabel($xaxstr)=1//, noLabel($yaxstr)=1
			//NameofPanel= "\\Z12"+stringfromlist(index-1,abcstr,";")+")"
			NameofPanel= "\\Z12Clst#"+num2istr(index)
			TextBox/C/F=0/B=1/N=$("name"+num2istr(index))/A=MC/X=((xhighFra-col_space/2)*100-50)/Y=(yhighFra*100-50) nameofpanel
			wvcount +=1
		endfor
	endfor

	
End

//This function plots cluster 1-x, plots only the individual thermograms on the same figure
//Each cluster ends up as a different panel in the same graph 
Function GraphClstTG0_OneGraph()

	string clusternumstr = "Clusternumsort"
	string ymatrixstr = "cmp_matrix_sn_10a_fil"
	string xwavestr = "temp_c_10a"
	variable ncols = 3
	
	prompt clusternumstr, "Select wave of cluster#", popup WaveList("*",";","")
	prompt ymatrixstr, "Select matrix of individual TGs", popup WaveList("*",";","")
	prompt xwavestr, "Select temperature wave", popup WaveList("*",";","")
	prompt ncols, "Columns of panels of graph"
	Doprompt "Select waves",clusternumstr, ymatrixstr, xwavestr, ncols
	
	if(V_flag)
		abort
	endif
	
	wave clusternum = $clusternumstr
	wave ymatrix = $ymatrixstr
	wave xwave = $xwavestr

	wavestats/q clusternum
	variable nclsts = V_max //number of clusters
	variable nions = dimsize(ymatrix,1) //number of individual ions
	variable nrows = ceil(nclsts/ncols) // number of rows
	
	colortab2wave dBZ21
	wave M_colors
	SetScale/I x,1,nclsts,M_colors
	variable red,green,blue
	
	variable col_offset
	variable col_len 
	variable col_space 
	variable col_multiplier = 2
	
	col_offset = 0.05
	//Let col_len=col_multiplier*col_space, then col_offset+ncols*col_len+(ncols-1)*col_space=1
	//so that
	col_space = (1-col_offset)/((1+col_multiplier)*ncols-1)
	col_len=col_multiplier*col_space
	
	variable row_len
	variable row_space = 0.2/nrows
	row_len = 1/nrows
	
	string yaxstr, xaxstr
	variable i, j,k
	variable wvcount = 0
	variable ioncount = 0
	variable index
	
	variable xlowFra,xhighFra,ylowFra,yhighFra
	string abcstr="a;b;c;d;e;f;g;h;i;j;k;l;m;n;o;p;q;r;s;t;u;v;w;x;y;z"
	string NameofPanel
	
	string graphname = "AvgThermogramsofClst"
	
	display/N=$graphname as graphname
	modifygraph width=500, height=500

	for(i=0;i<nrows;i+=1)
		for(j=0;j<ncols;j+=1)
			if(wvcount > nclsts-1)
				return 1
			endif
			index=j+1+ncols*i
			yaxstr = "y" + num2istr(index)
			xaxstr = "x" + num2istr(index)
			red=M_colors(index)[0]
			green=M_colors(index)[1]
			blue=M_colors(index)[2]
			for(k=0;k<nions;k+=1)
				if(clusternum[k]==index)
					appendtograph/l=$yaxstr/b=$xaxstr ymatrix[][k] vs xwave
					if(ioncount == 0)
						ModifyGraph rgb($(ymatrixstr))=(red,green,blue)
					else
						ModifyGraph rgb($(ymatrixstr+"#"+num2istr(ioncount)))=(red,green,blue)
					endif
					ioncount+=1
				endif
			endfor
	
			if(j==0)
				//ModifyGraph axisEnab($yaxstr)={1-row_len*(i+1)+row_space*(nrows-i),1-(row_len)*i-row_space*(i)},axisEnab($xaxstr)={col_offset+col_len*j,col_offset+col_len*(j+1)} // col_len*j+col_space*j,col_len2*(j+1)-col_space*(ncols-j-1)}
				ylowFra = 1-row_len*(i+1)+row_space //*(nrows-i)
				yhighFra = 1-(row_len)*i-row_space //*(i)
				xlowFra = col_offset+col_len*j
				xhighFra = col_offset+col_len*(j+1)
				ModifyGraph axisEnab($yaxstr)={ylowFra,yhighFra},axisEnab($xaxstr)={xlowFra,xhighFra}
			else
				ylowFra = 1-row_len*(i+1)+row_space //*(nrows-i)
				yhighFra = 1-(row_len)*i-row_space //*(i)
				xlowFra = col_offset+col_len*j+col_space*j
				xhighFra = col_offset+col_len*(j+1)+col_space*j
				ModifyGraph axisEnab($yaxstr)={ylowFra,yhighFra},axisEnab($xaxstr)={xlowFra,xhighFra}
			endif
			//print ylowFra,yhighFra,xlowFra,xhighFra
			SetAxis $yaxstr -0.1,1.3
			SetAxis $xaxstr 30,205			
			ModifyGraph freePos($yaxstr)={0,$xaxstr},freePos($xaxstr)={-0.1,$yaxstr}
			//ModifyGraph noLabel($xaxstr)=1//, noLabel($yaxstr)=1
			//NameofPanel= "\\Z12"+stringfromlist(index-1,abcstr,";")+")"
			NameofPanel= "\\Z12Clst#"+num2istr(index)
			TextBox/C/F=0/B=1/N=$("name"+num2istr(index))/A=MC/X=((xhighFra-col_space/2)*100-50)/Y=(yhighFra*100-50) nameofpanel
			wvcount +=1
		endfor
	endfor

	
End

Function GraphClstTG0_OneGraph_np(clusternum, ymatrix,xwave,ncols)

	wave clusternum
	wave ymatrix
	wave xwave
	variable ncols

	string clusternumstr=nameofwave(clusternum)
	string ymatrixstr=nameofwave(ymatrix)
	string xwavestr=nameofwave(xwave)

	wavestats/q clusternum
	variable nclsts = V_max //number of clusters
	variable nions = dimsize(ymatrix,1) //number of individual ions
	variable nrows = ceil(nclsts/ncols) // number of rows
	
	colortab2wave dBZ21
	wave M_colors
	SetScale/I x,1,nclsts,M_colors
	variable red,green,blue
	
	variable col_offset
	variable col_len 
	variable col_space 
	variable col_multiplier = 2
	
	col_offset = 0.05
	//Let col_len=col_multiplier*col_space, then col_offset+ncols*col_len+(ncols-1)*col_space=1
	//so that
	col_space = (1-col_offset)/((1+col_multiplier)*ncols-1)
	col_len=col_multiplier*col_space
	
	variable row_len
	variable row_space = 0.2/nrows
	row_len = 1/nrows
	
	string yaxstr, xaxstr
	variable i, j,k
	variable wvcount = 0
	variable ioncount = 0
	variable index
	
	variable xlowFra,xhighFra,ylowFra,yhighFra
	string abcstr="a;b;c;d;e;f;g;h;i;j;k;l;m;n;o;p;q;r;s;t;u;v;w;x;y;z"
	string NameofPanel
	
	string graphname = "AvgThermogramsofClst"
	
	display/N=$graphname as graphname
	modifygraph width=500, height=500

	for(i=0;i<nrows;i+=1)
		for(j=0;j<ncols;j+=1)
			if(wvcount > nclsts-1)
				return 1
			endif
			index=j+1+ncols*i
			yaxstr = "y" + num2istr(index)
			xaxstr = "x" + num2istr(index)
			red=M_colors(index)[0]
			green=M_colors(index)[1]
			blue=M_colors(index)[2]
			for(k=0;k<nions;k+=1)
				if(clusternum[k]==index)
					appendtograph/l=$yaxstr/b=$xaxstr ymatrix[][k] vs xwave
					if(ioncount == 0)
						ModifyGraph rgb($(ymatrixstr))=(red,green,blue)
					else
						ModifyGraph rgb($(ymatrixstr+"#"+num2istr(ioncount)))=(red,green,blue)
					endif
					ioncount+=1
				endif
			endfor
	
			if(j==0)
				//ModifyGraph axisEnab($yaxstr)={1-row_len*(i+1)+row_space*(nrows-i),1-(row_len)*i-row_space*(i)},axisEnab($xaxstr)={col_offset+col_len*j,col_offset+col_len*(j+1)} // col_len*j+col_space*j,col_len2*(j+1)-col_space*(ncols-j-1)}
				ylowFra = 1-row_len*(i+1)+row_space //*(nrows-i)
				yhighFra = 1-(row_len)*i-row_space //*(i)
				xlowFra = col_offset+col_len*j
				xhighFra = col_offset+col_len*(j+1)
				ModifyGraph axisEnab($yaxstr)={ylowFra,yhighFra},axisEnab($xaxstr)={xlowFra,xhighFra}
			else
				ylowFra = 1-row_len*(i+1)+row_space //*(nrows-i)
				yhighFra = 1-(row_len)*i-row_space //*(i)
				xlowFra = col_offset+col_len*j+col_space*j
				xhighFra = col_offset+col_len*(j+1)+col_space*j
				ModifyGraph axisEnab($yaxstr)={ylowFra,yhighFra},axisEnab($xaxstr)={xlowFra,xhighFra}
			endif
			//print ylowFra,yhighFra,xlowFra,xhighFra
			SetAxis $yaxstr -0.1,1.3
			SetAxis $xaxstr 30,205			
			ModifyGraph freePos($yaxstr)={0,$xaxstr},freePos($xaxstr)={-0.1,$yaxstr}
			//ModifyGraph noLabel($xaxstr)=1//, noLabel($yaxstr)=1
			//NameofPanel= "\\Z12"+stringfromlist(index-1,abcstr,";")+")"
			NameofPanel= "\\Z12Clst#"+num2istr(index)
			TextBox/C/F=0/B=1/N=$("name"+num2istr(index))/A=MC/X=((xhighFra-col_space/2)*100-50)/Y=(yhighFra*100-50) nameofpanel
			wvcount +=1
		endfor
	endfor

	
End


///////////////////////////////////////////////////////////////////////////////////////////////////
//This function graphs the spectrum of each cluster in a stacked way
Function GraphClstSpectrum_stacked()
	
	string Spectramatstr = "ClusterSpectra_normsort"
	string masswavestr = "Molarmass_cmp_fil"
	
	prompt Spectramatstr, "Select matrix of spectra", popup WaveList("*",";","")
	prompt masswavestr, "Select m/z wave", popup WaveList("*",";","")
	Doprompt "Select waves",Spectramatstr, masswavestr
	
	if(V_flag)
		abort
	endif
	
	wave Spectramat = $Spectramatstr
	wave masswave = $masswavestr
	
	variable nclsts,i
	string name

	nclsts=dimsize(Spectramat,0)
	
	colortab2wave dBZ21
	wave M_colors
	SetScale/I x,1,nclsts,M_colors
	variable red,green,blue
   
	string graphname = "StackedSpectrum"
	
	display/N=$graphname as graphname
	modifygraph width=300, height=350
	
	for(i=0;i<nclsts;i+=1)
	
		red=M_colors(i+1)[0]
		green=M_colors(i+1)[1]
		blue=M_colors(i+1)[2]
		
		name="Clst"+num2istr(i+1)
		appendtograph/l=$name Spectramat[i][] vs masswave
		
		if(i==0)
			ModifyGraph rgb($(nameofwave(Spectramat)))=(red,green,blue)
		else
			ModifyGraph rgb($(nameofwave(Spectramat)+"#"+num2istr(i)))=(red,green,blue)
		endif
		
	endfor
   
	modifygraph mode=1,lsize=1.5
	SetAxis bottom 10,*
	stackallaxes_cdc("",0,0)
   
   	string nameofpanel
   	
   	variable yspace=1/nclsts
   	
  	for(i=0;i<nclsts;i+=1)
  		name="Clst"+num2str(i+1)
  		
		modifygraph zero($name)=1,freepos($name)={10,bottom}
		Setaxis $name 0,*
		ModifyGraph mirror($name)=1,standoff($name)=0
		
		NameofPanel= "\\Z12Clst#"+num2istr(i+1)
		TextBox/C/F=0/B=1/N=$("name"+num2istr(i+1))/A=MC/X=(-35)/Y=(((i+1)*yspace)*100-53) nameofpanel
	endfor
   
   	TextBox/C/N=$("Lbl_left")/O=90/F=0/B=1/A=LC/X=(-20)/Y=0 "Normalized mass"
	Label bottom "m/z"
	//Label left "Absolute Signal"
   
  	ModifyGraph mirror(bottom)=1,standoff(bottom)=0
  	//ModifyGraph axisEnab(bottom)={0.2,0.95}
   
End


Function GraphClstSpectrum_stacked_np(Spectramat,masswave)
	wave Spectramat //2D wave
	wave masswave
	
	variable nclsts,i
	string name

	nclsts=dimsize(Spectramat,0)
	
	colortab2wave dBZ21
	wave M_colors
	SetScale/I x,1,nclsts,M_colors
	variable red,green,blue
   
	string graphname = "StackedSpectrum"
	
	display/N=$graphname as graphname
	modifygraph width=300, height=350
	
	for(i=0;i<nclsts;i+=1)
	
		red=M_colors(i+1)[0]
		green=M_colors(i+1)[1]
		blue=M_colors(i+1)[2]
		
		name="Clst"+num2istr(i+1)
		appendtograph/l=$name Spectramat[i][] vs masswave
		
		if(i==0)
			ModifyGraph rgb($(nameofwave(Spectramat)))=(red,green,blue)
		else
			ModifyGraph rgb($(nameofwave(Spectramat)+"#"+num2istr(i)))=(red,green,blue)
		endif
		
	endfor
   
	modifygraph mode=1,lsize=1.5
	SetAxis bottom 10,*
	stackallaxes_cdc("",0,0)
   
   	string nameofpanel
   	
   	variable yspace=1/nclsts
   	
  	for(i=0;i<nclsts;i+=1)
  		name="Clst"+num2str(i+1)
  		
		modifygraph zero($name)=1,freepos($name)={10,bottom}
		Setaxis $name 0,*
		ModifyGraph mirror($name)=1,standoff($name)=0
		
		NameofPanel= "\\Z12Clst#"+num2istr(i+1)
		TextBox/C/F=0/B=1/N=$("name"+num2istr(i+1))/A=MC/X=(-35)/Y=(((i+1)*yspace)*100-53) nameofpanel
	endfor
   
   	TextBox/C/N=$("Lbl_left")/O=90/F=0/B=1/A=LC/X=(-15)/Y=0 "Normalized mass"
	Label bottom "m/z"
	//Label left "Absolute Signal"
   
  	ModifyGraph mirror(bottom)=1,standoff(bottom)=0
  	//ModifyGraph axisEnab(bottom)={0.2,0.95}
   
End


// Usage Notes for StackAllAxes
// StackAllAxes looks for all left axes and stacks them
// on graphName
// lo_twink and hi_twink are the number of hundredth to inwardly twink the middle axes by
Function StackAllAxes_cdc(graphName, lo_twink, hi_twink)
	String graphName
	Variable lo_twink, hi_twink
	
	String axes_list = AxisList( graphName ), this_ax, this_ax_type
	axes_list = RemoveFromList(  "top", axes_list)
	axes_list = RemoveFromList( "bottom", axes_list)
	Variable num_axes = ItemsInList( axes_list )
	
	if( num_axes == 0 )
		return -1
	endif
	Variable kdex = 0
	do
		this_ax = StringFromList( kdex, axes_list )
		this_ax_type = StringByKey( "AXTYPE", AxisInfo( graphName, this_ax ) )
		if( (cmpstr( lowerstr( this_ax_type), "top" ) == 0)      %|     (cmpstr( lowerstr( this_ax_type), "bottom" ) == 0) )
			axes_list = RemoveFromList( this_ax, axes_list )
		endif
		kdex += 1
	while( kdex < num_axes )
		 num_axes = ItemsInList( axes_list )
	Variable idex = 0, interval = 1/num_axes, low_bound, high_bound
	
	
	do
		this_ax = StringFromList( idex, axes_list )
		low_bound = idex * interval
		high_bound = (idex + 1) * interval
		if( (idex >= 0) %& (idex < num_axes + 1 ) )
	//		low_bound += lo_twink/100;	high_bound -= hi_twink/100 // changed by CD Cappa 07/30/10
			low_bound = low_bound; 		high_bound -= hi_twink/100
		endif
		if( low_bound < 0 )
			low_bound = 0
		endif
		if( high_bound > 1 )
			high_bound = 1 
		endif
		ModifyGraph axisEnab($this_ax)={low_bound,high_bound}
		idex += 1
	while( idex < num_axes )
	return 1
	
End

//This function graphs the stacked sum thermograms of clusters, including the bulk thermogram	
Function GraphClstTG_stacked()

	string ymatrixstr = "SumTG_clstsort"
	string ywavestr = "Bulk_TG_rf"
	string xwavestr = "time_s"

	prompt ymatrixstr, "Select matrix of individual TGs", popup WaveList("*",";","")
	prompt ywavestr, "Select bulk TG wave", popup WaveList("*",";","")
	prompt xwavestr, "Select desorption time wave", popup WaveList("*",";","")
	Doprompt "Select waves",ymatrixstr, ywavestr, xwavestr
	
	if(V_flag)
		abort
	endif

	wave ymatrix = $ymatrixstr //Matrix of sum TG of clusters
	wave ywave = $ywavestr //Bulk TG
	wave xwave = $xwavestr	 //Usually time wave
	
	string graphname = "StackedThermogram"
	
	GraphMatrix_Col_np(ymatrix, xwave, graphname)
	
	ModifyGraph mode=7,toMode=2
	ModifyGraph hbFill=2
	appendtograph ywave vs xwave
	ReorderTraces $ymatrixstr,{$ywavestr}
	ModifyGraph lsize($ywavestr)=1.5,rgb($ywavestr)=(0,0,0)
	ModifyGraph mirror=1,standoff=0
	Label left "Absolute mass"
	Label bottom "Desorption time, s"
	Legend/C/N=text0/J/A=RT "\\s("+ywavestr+") Bulk Thermogram"
	
	modifygraph width=250, height=150 

End

Function GraphClstTG_stacked_np(ymatrix,ywave,xwave)

	wave ymatrix //Matrix of sum TG of clusters
	wave ywave //Bulk TG
	wave xwave	 //Usually time wave
	
	string graphname = "StackedThermogram"
	
	GraphMatrix_Col_np(ymatrix, xwave,graphname)
	
	string ymatrixstr = nameofwave(ymatrix)
	string ywavestr = nameofwave(ywave)
	string xwavestr = nameofwave(xwave)
	
	ModifyGraph mode=7,toMode=2
	ModifyGraph hbFill=2
	appendtograph ywave vs xwave
	ReorderTraces $ymatrixstr,{$ywavestr}
	ModifyGraph lsize($ywavestr)=1.5,rgb($ywavestr)=(0,0,0)
	ModifyGraph mirror=1,standoff=0
	Label left "Absolute mass"
	Label bottom "Desorption time, s"
	Legend/C/N=text0/J/A=RT "\\s("+ywavestr+") Bulk Thermogram"
	
	modifygraph width=250, height=150 

End
	
//Graph each column vs. xwave for all the columns
Function GraphMatrix_Col_np(ymatrix, xwave, graphname)
	wave ymatrix, xwave
	string graphname
	
	variable ncols=dimsize(ymatrix,1)
	variable nrows=dimsize(ymatrix,0)
	
	if(nrows!=numpnts(xwave))
		print "The number of rows of the matrix doesn't match the number points of the xwave"
		abort
	endif
	
	colortab2wave dBZ21
	wave M_colors
	SetScale/I x,1,ncols,M_colors
	variable red,green,blue
	
	make/o/d/n=(ncols,3) ClstColors
	
	variable i
	for(i=0 ; i<ncols ; i+=1)
	
		red=M_colors(i+1)[0]
		green=M_colors(i+1)[1]
		blue=M_colors(i+1)[2]
		
		ClstColors[i][0] = red
		ClstColors[i][1] = green
		ClstColors[i][2] = blue		
		
		if(i==0)
			display/N=$graphname as graphname
			appendtograph ymatrix[][i] vs xwave
			ModifyGraph rgb($(nameofwave(ymatrix)))=(red,green,blue)
		else
			appendtograph ymatrix[][i] vs xwave
			ModifyGraph rgb($(nameofwave(ymatrix)+"#"+num2istr(i)))=(red,green,blue)
		endif
	endfor
End

//This function graphs the relationship of noise level and mass fractional contribution
Function GraphNoiseMass_np(Noisewave, masswave,filterwave)

	wave noisewave, masswave, filterwave
	
	string noisewavestr = nameofwave(noisewave)
	string masswavestr = nameofwave(masswave)
	string filterwavestr = nameofwave(filterwave)
	string graphname = "Noise_v_Mass"
	
	display/N=$graphname as graphname
	modifygraph width=250,height=200
	
	appendtograph noisewave vs masswave
	
	ModifyGraph mode=3
	ModifyGraph log=1
	Label left "Noise level"
	Label bottom "Mass fractional contribution"
	ModifyGraph zColor($noisewavestr)={$filterwavestr,*,*,CyanMagenta,0}
	TextBox/C/N=text0/A=RT "Cyan -- Screened out\rPink -- Kept for clustering"

End

Function GraphMassThreshold_np(cmp_mass_norm_sort, cmp_mass_norm_int)

	wave cmp_mass_norm_sort, cmp_mass_norm_int
	
	display cmp_mass_norm_sort vs cmp_mass_norm_int
	modifygraph width=250,height=200
	
	ModifyGraph lowTrip=0.01
	SetAxis bottom 0,1.01
	//SetAxis left -0.01, *
	ModifyGraph mirror=1,fSize=12,standoff=0
	ModifyGraph mode=3,marker=19
	Label left "Individual mass contribution"
	Label bottom "Accumulative mass contribution"
	
End
	
//This function is usually used after executing function "FindOptEps_NSSC"
//The Cluster# of ions are stored in the columns of a 2D matrix, "Clusternummatrix"
//Choose a specific epsilon to plot the  clustered thermograms
Function GraphTG0_fromMat()
	
	string clusternummatrixstr = "Clusternummatrix"
	string epswavestr = "epsilon"
	string ymatrixstr = "cmp_matrix_sn_10a_fil"
	string xwavestr = "temp_c_10a"
	variable epsindx = 10
	
	prompt clusternummatrixstr, "Select matrix of cluster#", popup WaveList("*",";","")
	prompt epswavestr, "Select epsilon wave", popup WaveList("*",";","")
	prompt ymatrixstr, "Select matrix of individual TGs", popup WaveList("*",";","")
	prompt xwavestr, "Select time or T wave", popup WaveList("*",";","")
	prompt epsindx, "Enter the index # of epsilon"
	Doprompt "Select waves",clusternummatrixstr,epswavestr,ymatrixstr, xwavestr, epsindx
	
	if(V_flag)
		abort
	endif
	
	wave clusternummatrix = $clusternummatrixstr
	wave epswave = $epswavestr
	wave ymatrix = $ymatrixstr
	wave xwave = $xwavestr
	
	variable nions=dimsize(ymatrix,1)
	
	make/o/d/n=(nions)/Free tempClusternum
	
	tempclusternum=clusternummatrix[p][epsindx] 
	
	wavestats/q tempclusternum
	variable nclusters=V_max
	variable ncols=ceil(nclusters^0.5)
	
	GraphClstTG0_OneGraph_np(tempclusternum, ymatrix,xwave,ncols)
	
	TextBox/C/F=0/N=text0/A=MC/Y=50 "Eps = "+num2str(epswave[epsindx])
End


Function GraphTG0_fromMat_np(Clusternummatrix,epsindx,epswave,ymatrix,xwave)
	wave clusternummatrix
	variable epsindx
	wave epswave
	wave ymatrix
	wave xwave
	
	variable nions=dimsize(ymatrix,1)
	
	make/o/d/n=(nions)/Free tempClusternum
	
	tempclusternum=clusternummatrix[p][epsindx] 
	
	wavestats/q tempclusternum
	variable nclusters=V_max
	variable ncols=ceil(nclusters^0.5)
	
	GraphClstTG0_OneGraph_np(tempclusternum, ymatrix,xwave,ncols)
	
	TextBox/C/F=0/N=text0/A=MC/Y=50 "Eps = "+num2str(epswave[epsindx])
End


Function GraphTG0_fromMat_range_np(Clusternummatrix,starteps,endeps,epswave,ymatrix,xwave)
	
	wave clusternummatrix
	variable starteps, endeps
	wave epswave
	wave ymatrix
	wave xwave
	
	variable i
	
	for(i=starteps;i<=endeps;i+=1)
		GraphTG0_fromMat_np(Clusternummatrix,i,epswave,ymatrix,xwave)
	endfor
	
End

//This function graphs the size (number of ions in a cluster) of clusters as bar chart
Function GraphClstSize()
	
	string clustersizestr = "Clustersizesort"
	
	prompt clustersizestr, "Select the wave of cluster size",  popup WaveList("*",";","")
	Doprompt "Select waves", clustersizestr

	if(V_flag)
		abort
	endif
	
	wave clustersizewave = $clustersizestr

	variable nclusters = numpnts(clustersizewave)
	wavestats/q clustersizewave
	variable maxsize = V_max
	
	make/o/d/n=(nclusters) Clusterindx
	Clusterindx = x
	Note/k Clusterindx "This wave stores the index# of clusters from 1 to "+num2istr(nclusters-1)+" ; point 0 indicates unclustered"
	
	string Graphname = "ClusterSize"
	
	display/N=$graphname as graphname
	modifygraph width=350, height=150
	appendtograph clustersizewave vs clusterindx
	
	ModifyGraph mode=8,lsize=10
	ModifyGraph textMarker($clustersizestr)={$clustersizestr,"default",0,0,5,0.00,10.00}
	ModifyGraph zColor($clustersizestr)={Clusterindx,1,*,dBZ21,0}
	ModifyGraph zColorMin($Clustersizestr)=(0,0,0)
	ModifyGraph mirror=1,standoff(bottom)=0
	ModifyGraph nticks(bottom)=nclusters
	Label left "Number of thermograms"
	Label bottom "Cluster # (0 = unclustered)"
	ModifyGraph nticks(bottom)=10
	SetAxis left 0,(maxsize+5)
	
End


Function GraphClstSize_np(clustersizewave)

	wave clustersizewave

	variable nclusters = numpnts(clustersizewave)
	wavestats/q clustersizewave
	variable maxsize = V_max
	
	string clustersizestr = nameofwave(clustersizewave)
	
	make/o/d/n=(nclusters) Clusterindx
	Clusterindx = x
	Note/k Clusterindx "This wave stores the index# of clusters from 1 to "+num2istr(nclusters-1)+" ; point 0 indicates unclustered"
	
	string Graphname = "ClusterSize"
	
	display/N=$graphname as graphname
	modifygraph width=350, height=150
	appendtograph clustersizewave vs clusterindx
	
	ModifyGraph mode=8,lsize=10
	ModifyGraph textMarker($clustersizestr)={$clustersizestr,"default",0,0,5,0.00,10.00}
	ModifyGraph zColor($clustersizestr)={Clusterindx,1,*,dBZ21,0}
	ModifyGraph zColorMin($Clustersizestr)=(0,0,0)
	ModifyGraph mirror=1,standoff(bottom)=0
	ModifyGraph nticks(bottom)=nclusters
	Label left "Number of thermograms"
	Label bottom "Cluster # (0 = unclustered)"
	ModifyGraph nticks(bottom)=10
	SetAxis left 0,(maxsize+5)
	
End

//This function graphs the mass fractional contribution of clusters as bar chart
Function GraphClstMass()

	string clustermassstr = "Clustermasssort"
	
	prompt clustermassstr, "Select the wave of cluster mass fraction",  popup WaveList("*",";","")
	Doprompt "Select waves", clustermassstr

	if(V_flag)
		abort
	endif
	
	wave clustermasswave = $clustermassstr
	
	variable nclusters = numpnts(clustermasswave)

	make/o/d/n=(nclusters+1) Clusterindx_ext
	Clusterindx_ext = x-1
	
	Note/k Clusterindx_ext "This wave stores the index# of clusters from 2 to "+num2istr(nclusters)+" ; point 0 stores filtered and point 1 stores unclustered"
	
	make/o/d/n=(nclusters) $(clustermassstr+"_ext")
	wave clustermassext = $(clustermassstr+"_ext")
	
	Note/k Clustermassext "This wave stores the percentage of mass contribution of filtered ions (point 0); unclustered ions (point 1) and clustered ions (point 2 - "+num2istr(nclusters)
	
	string clustermassextstr = clustermassstr+"_ext" 
	clustermassext = clustermasswave
	InsertPoints 0,1, Clustermassext
	wavestats/q clustermassext
	Clustermassext[0]=1-V_sum
	Clustermassext*=100
	
	wavestats/q clustermassext
	variable maxmass = V_max
	
	string Graphname = "ClusterMassFra"
	
	display/N=$graphname as graphname
	modifygraph width=400, height=150
	
	appendtograph clustermassext vs clusterindx_ext
	
	ModifyGraph mode=8,lsize=10
	ModifyGraph zColor($Clustermassextstr)={Clusterindx_ext,1,*,dBZ21,0}
	ModifyGraph zColorMin($Clustermassextstr)=(0,0,0)
	ModifyGraph mirror=1,standoff(bottom)=0
	ModifyGraph nticks(bottom)=nclusters+1
	Label left "Mass contribution, %"
	Label bottom "Cluster # (0 = unclustered, -1 = filtered out)"
	ModifyGraph textMarker($Clustermassextstr)={$Clustermassextstr,"default",0,0,5,0.00,5.00}
	SetAxis left 0, maxmass+5
	ModifyGraph msize=3

End


Function GraphClstMass_np(Clustermasswave)
	wave clustermasswave

	string clustermassstr = nameofwave(clustermasswave)
	
	variable nclusters = numpnts(clustermasswave)

	make/o/d/n=(nclusters+1) Clusterindx_ext
	Clusterindx_ext = x-1
	
	Note/k Clusterindx_ext "This wave stores the index# of clusters from 2 to "+num2istr(nclusters)+" ; point 0 stores filtered and point 1 stores unclustered"
	
	make/o/d/n=(nclusters) $(clustermassstr+"_ext")
	wave clustermassext = $(clustermassstr+"_ext")
	
	Note/k Clustermassext "This wave stores the percentage of mass contribution of filtered ions (point 0); unclustered ions (point 1) and clustered ions (point 2 - "+num2istr(nclusters)
	
	string clustermassextstr = clustermassstr+"_ext" 
	clustermassext = clustermasswave
	InsertPoints 0,1, Clustermassext
	wavestats/q clustermassext
	Clustermassext[0]=1-V_sum
	Clustermassext*=100
	
	wavestats/q clustermassext
	variable maxmass = V_max
	
	string Graphname = "ClusterMassFra"
	
	display/N=$graphname as graphname
	modifygraph width=400, height=150
	
	appendtograph clustermassext vs clusterindx_ext
	
	ModifyGraph mode=8,lsize=10
	ModifyGraph zColor($Clustermassextstr)={Clusterindx_ext,1,*,dBZ21,0}
	ModifyGraph zColorMin($Clustermassextstr)=(0,0,0)
	ModifyGraph mirror=1,standoff(bottom)=0
	ModifyGraph nticks(bottom)=nclusters+1

	Label left "Mass contribution, %"
	Label bottom "Cluster # (0 = unclustered, -1 = filtered out)"
	ModifyGraph textMarker($Clustermassextstr)={$Clustermassextstr,"default",0,0,5,0.00,5.00}
	SetAxis left 0, maxmass+5
	ModifyGraph msize=3

End

//This function converts location (index#) of a wave to the Temperature according to a T wave
Function Loc2T_np(locwave,Twave,newwavestr)
	wave locwave, Twave
	string newwavestr
	
	variable i
	
	variable nrows=numpnts(locwave)
	
	make/o/d/n=(nrows) $newwavestr
	wave newwave=$newwavestr
	
	for(i=0;i<nrows;i+=1)
		newwave[i]=Twave[locwave[i]]
	endfor
	
	Note/k newwave "This wave stores the values from wave "+nameofwave(Twave)+" according to index in wave "+nameofwave(locwave)
	
End


Function GraphBasics()

	wave cmp_mass_rf
	wave molarmass_cmp_rf
	
	display cmp_mass_rf vs molarmass_cmp_rf
	ModifyGraph width=350, height=150
	ModifyGraph mode=1
	ModifyGraph mirror=1,fSize=12,standoff=0
	Label left "Mass (sum of mass over \rentire desorption time)"
	Label bottom "m/z (without iodine)"
	
	wave Bulk_TG_rf
	wave temp_c, time_s
	
	display bulk_TG_rf vs temp_c
	ModifyGraph width=250, height=200
	Label left "Mass (sum of mass of all the ions)"
	Label bottom "Desorption temperature, ¡C"
	ModifyGraph lsize=1.5
	ModifyGraph mirror=1,fSize=12,standoff=0
	Setaxis left 0,*
	
	display bulk_TG_rf vs time_s
	appendtograph/r temp_c vs time_s
	ModifyGraph width=250, height=200
	Label left "Mass (sum of mass of all the ions)"
	Label bottom "Desorption time, s"
	ModifyGraph lsize=1.5
	ModifyGraph rgb(temp_C)=(32768,40777,65535)
	ModifyGraph lstyle(temp_C)=3
	Label right "Desorption temperature, ¡C"
	ModifyGraph axRGB(right)=(1,4,52428),tlblRGB(right)=(1,4,52428)
	ModifyGraph alblRGB(right)=(1,4,52428)
	ModifyGraph mirror(bottom)=1,fSize=12,standoff=0
End

Function GraphmultiTG_np(ywvstr,xwvstr)
	string ywvstr  //By default refers to a 2D matrix with columns containing TGs of clusters
	string xwvstr
	
	string ywvlist = wavelist(ywvstr+"*",";","")
	string xwvlist = wavelist(xwvstr+"*",";","")
	variable nwvs = itemsinlist(ywvlist)
	
	colortab2wave BlueGreenOrange
	wave M_colors
	SetScale/I x,0,nwvs-1,M_colors
	variable red,green,blue
	
	string ywvlistsort = sortlist(ywvlist,";",16)
	string xwvlistsort = sortlist(xwvlist,";",16)
		
	wave tempwv = $(stringfromlist(0,ywvlistsort))
	variable nclsts = dimsize(tempwv,1)
	
	variable ncols = ceil((nclsts)^0.5)  //number of cols of panels
	variable nrows = ceil(nclsts/ncols) //number of rows of panels
	
	variable col_offset
	variable col_len 
	variable col_space 
	variable col_multiplier = 2
	
	col_offset = 0.05
	col_space = (1-col_offset)/((1+col_multiplier)*ncols-1)
	col_len=col_multiplier*col_space
	
	variable row_len
	variable row_space = 0.2/nrows
	row_len = 1/nrows
	
	string yaxstr, xaxstr
	variable i,j,k
	variable wvcount = 0
	variable index
	
	variable xlowFra,xhighFra,ylowFra,yhighFra
	string NameofPanel
	
	string graphname = ywvstr+"_v_"+xwvstr
	
	display/N=$graphname as graphname
	modifygraph width=500, height=500
	
	string tempywvstr, tempxwvstr

	for(i=0;i<nrows;i+=1)
		for(j=0;j<ncols;j+=1)
			if(wvcount > nclsts-1)
				return 1
			endif
			index=j+1+ncols*i
			yaxstr = "y" + num2istr(index)
			xaxstr = "x" + num2istr(index)
			for(k=0;k<nwvs;k+=1)
				tempywvstr = stringfromlist(k,ywvlistsort)
				tempxwvstr = stringfromlist(k,xwvlistsort)
				wave tempywv = $tempywvstr
				wave tempxwv = $tempxwvstr
				red=M_colors(k)[0]
				green=M_colors(k)[1]
				blue=M_colors(k)[2]
				appendtograph/l=$yaxstr/b=$xaxstr tempywv[][wvcount] vs tempxwv
				if(wvcount == 0)
					ModifyGraph rgb($(tempywvstr))=(red,green,blue)
				else
					ModifyGraph rgb($(tempywvstr+"#"+num2istr(wvcount)))=(red,green,blue)
				endif
			endfor
			
			if(j==0)
				//ModifyGraph axisEnab($yaxstr)={1-row_len*(i+1)+row_space*(nrows-i),1-(row_len)*i-row_space*(i)},axisEnab($xaxstr)={col_offset+col_len*j,col_offset+col_len*(j+1)} // col_len*j+col_space*j,col_len2*(j+1)-col_space*(ncols-j-1)}
				ylowFra = 1-row_len*(i+1)+row_space //*(nrows-i)
				yhighFra = 1-(row_len)*i-row_space //*(i)
				xlowFra = col_offset+col_len*j
				xhighFra = col_offset+col_len*(j+1)
				ModifyGraph axisEnab($yaxstr)={ylowFra,yhighFra},axisEnab($xaxstr)={xlowFra,xhighFra}
			else
				ylowFra = 1-row_len*(i+1)+row_space //*(nrows-i)
				yhighFra = 1-(row_len)*i-row_space //*(i)
				xlowFra = col_offset+col_len*j+col_space*j
				xhighFra = col_offset+col_len*(j+1)+col_space*j
				ModifyGraph axisEnab($yaxstr)={ylowFra,yhighFra},axisEnab($xaxstr)={xlowFra,xhighFra}
			endif
			//print ylowFra,yhighFra,xlowFra,xhighFra
			//SetAxis $yaxstr -0.1,1.3
			//SetAxis $xaxstr 30,205			
			ModifyGraph freePos($yaxstr)={0,$xaxstr},freePos($xaxstr)={0,$yaxstr}
			ModifyGraph lsize=1.5

			NameofPanel= "\\Z12Clst#"+num2istr(index)
			TextBox/C/F=0/B=1/N=$("name"+num2istr(index))/A=MC/X=((xhighFra-col_space/2)*100-50)/Y=(yhighFra*100-50) nameofpanel
			wvcount +=1
		endfor
	endfor
	
	
End


Function SetAxisRange_panels(axstr,lowv,highv,nclsts)
	string axstr
	variable lowv,highv
	variable nclsts
	
	variable i
	string yaxstr
	
	for(i=0;i<nclsts;i+=1)
		yaxstr = axstr + num2istr(i+1)

		SetAxis $yaxstr lowv,highv
	endfor
	
End


Function Calculate_RinterClst_indx(indxwave,ymatrix,eps)
	wave indxwave
	wave ymatrix
	variable eps
	
	variable nrows = dimsize(ymatrix,0)
	variable nclstN = numpnts(indxwave)
	
	variable j
	
	make/o/d/n=(nrows,nclstN)/Free clusterseed_matrix
			
	for(j=0;j<nclstN;j+=1)
		clusterseed_matrix[][j]=ymatrix[p][indxwave[j]]
	endfor
	
	CalculateED_np(clusterseed_matrix,"ED_clsd")
	wave waveclusterSD=$("ED_clsd")
	wavestats/q waveclusterSD
	
	variable R_interClst
	R_interClst=V_sum/(nClstN*(nClstN-1))/eps
	
	killwaves waveclusterSD
	
	return R_interClst	
End

Function Calculate_RinterClst_TG(TGmatrix,eps)
	wave TGmatrix
	variable eps
	
	variable nrows = dimsize(TGmatrix,0)
	variable nclstN = dimsize(TGmatrix,1)
	
	variable j
	
	CalculateED_np(TGmatrix,"ED_clsd")
	wave waveclusterSD=$("ED_clsd")
	wavestats/q waveclusterSD
	
	variable R_interClst
	R_interClst=V_sum/(nClstN*(nClstN-1))/eps
	
	killwaves waveclusterSD
	
	return R_interClst	
End


Function CalculateTm50(Clustermatrix,xwave,massfra,newwavestr)
	wave clustermatrix
	wave xwave
	variable massfra
	string newwavestr
	
	variable nrows = dimsize(clustermatrix,0)
	variable ncols = dimsize(clustermatrix,1)
	make/o/d/n=(ncols) $(newwavestr)
	wave newwave = $(newwavestr)
	make/o/d/n=(nrows)/Free tempwave1
	
	make/o/d/n=(nrows) $(nameofwave(xwave)+"_im")
	wave xwave_im=$(nameofwave(xwave)+"_im")
	xwave_im=xwave
	Insertpoints nrows,1,xwave_im
	xwave_im[nrows]=xwave_im[nrows-1]+(xwave[nrows-1]-xwave[nrows-2])
	
	variable i,j
	variable fitstartpnt
	
	for(i=0;i<ncols;i+=1)
		tempwave1=clustermatrix[p][i]
		wavestats/q tempwave1
		Integrate tempwave1/X=xwave_im/D=temp_INT
		wave temp_INT
		for(j=0;j<nrows-1;j+=1)
			if(temp_INT[j]<=massfra*temp_INT[nrows-1] && temp_INT[j+1]>=massfra*temp_INT[nrows-1])
				newwave[i]=xwave[j]
			endif
		endfor
	endfor
	
	killwaves xwave_im, temp_int
End