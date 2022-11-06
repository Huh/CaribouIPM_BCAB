Barry Nobert
16 July 2021

Caribou_Demographic_Vital_Rates_9July2021 - 2021 cairbou demographic data to be put in the cairbou database
as a table. This way the data is secure and the latest/update is apparent. See column definitions below:

	Reporting_Pop -- population that demographics are estimated for

	Range_ID -- unique identifier of the Range within the caribou db

	Herd_ID -- unique identifier of the Herd within the caribou db (herd is nested within range)

	Year_END -- caribou year end (April of each year).          
 
	Recruitment_Mean --  mean of annual female-only recruitment rate per adult female

	Recruitment_SD --  standard deviation of annual female-only recruitment rate per adult female

	Recruitment_LCL95  --  lower 95% confidence interval of annual female-only recruitment rate per adult female

	Recruitment_UCL95 --  upper 95% confidence interval of annual female-only recruitment rate per adult female

	Survival_Mean       --    mean of adult female annual survival rate

	Survival_SD       --    standard deviation of adult female annual survival rate

	Survival_LCL95      --    lower 95% confidence interval of adult female annual survival rate
 
	Survival_UCL95       --  upper 95% confidence interval of adult female annual survival rate
 
	Lambda_Mean       --  mean of adult female population growth rate

	Lambda_SD       --  standard deviation of adult female population growth rate

	Lambda_LCL95      --  lower 95% confidence interval of adult female population growth rate

	Lambda_UCL95      --   upper 95% confidence interval of adult female population growth rate

	DemoMethod         --  method used to estimate demographics (Eaker = Shiny R app, Hervieux = PopTools in excel). 

	LastUpdated        -- date that demographics were last updated

	UpdatedBy          --  initials of staff that did the demographic update

	Treatment         --  whether wolf control happened or not

	KM_mean           --   mean emperical Kaplan–Meier estimate of adult female survival

	KM_se             --  standard error emperical Kaplan–Meier estimate of adult female survival 

	CollarsStartofYear --  defacto sample size for the KM estimates (# of collars at the start of the caribou year in May)

	CSTotGroups         --  calf survey - number of groups observed

	CSTotCollars       --   calf survey - number of collared animals observed (includes dead collars) 

	CSTotCaribou       --  calf survey - number of caribou observed

	CSTotCows         --   calf survey - number of adult/yearling females observed

	CSTotBulls        -- calf survey - number of adult/yearling bulls observed

	CSTotCalves       --  calf survey - number of calves (either sex) observed

	Comments	-- additional clarification of the reporting populations

