# Stressing step
# Mechanical layering
# Disp. BC: Top Disp_y = 0; left Disp_x = 0; right Disp_x -> function
# Constrain - Slider: Bot Disp_y
# Temp. BC: Top = function depth / constant, Bot = function depth / constant
# IC: Temp- function depth;
# Executioner with "Returnmap" - timestepper
[GlobalParams]
  time_factor = 1.e-3
[]

[Mesh]
  type = FileMesh
  file = Layered_Fault_Real_fine_quad.msh
[]

[Variables]
  active = 'temp disp_y disp_x pore_pressure'
  [./disp_x]
    order = FIRST
    family = LAGRANGE
  [../]
  [./disp_y]
    order = FIRST
    family = LAGRANGE
  [../]
  [./temp]
  [../]
  [./pore_pressure]
  [../]
[]

[Materials]
[./Fault_nomech]
   type = RedbackMaterial
   block = '5'
   gr = 1e-4
   ar = 10
   temperature = temp
   pore_pres = pore_pressure
   phi0 = 0.1
   Kc = 10
   ref_lewis_nb = 1
   ar_F = 16
   disp_x = disp_x
   ar_R = 10
   disp_y = disp_y
   eta1 = 1e3
   mu = 1
   Aphi = 1
   is_mechanics_on = true
   # total_porosity = total_porosity
   solid_thermal_expansion = 1e-3
   #gravity = '0 -9.81 0'
   solid_density = 2.5
   solid_compressibility = 1.1111111111
  [../]
  [./Fault_mech]
    type = RedbackMechMaterialJ2
    block = '5'
    youngs_modulus = 5
    exponent = 1
    disp_y = disp_y
    disp_x = disp_x
    yield_stress = '0. 0.1 1. 0.1'
    temperature = temp
    pore_pres = pore_pressure
    poisson_ratio = 0.25
    ref_pe_rate = 1
    # outputs = all
  [../]
[./Quartaer_nomech]
    type = RedbackMaterial
    block = '1'
    gr = 1e-5
    ar = 12
    temperature = temp
    pore_pres = pore_pressure
    phi0 = 0.1
    Kc = 10
    ref_lewis_nb = 1
    ar_F = 16
    disp_x = disp_x
    ar_R = 10
    disp_y = disp_y
    eta1 = 1e3
    mu = 1
    Aphi = 1
    is_mechanics_on = true
    # total_porosity = total_porosity
    solid_thermal_expansion = 1e-3
    #gravity = '0 -9.81 0'
    solid_density = 2.65
    solid_compressibility = 1.1111111111
  [../]
  [./Quartaer_mech]
    type = RedbackMechMaterialJ2
    block = '1'
    youngs_modulus = 15
    exponent = 5
    disp_y = disp_y
    disp_x = disp_x
    yield_stress = '0. 1 1. 1'
    temperature = temp
    pore_pres = pore_pressure
    poisson_ratio = 0.3
    ref_pe_rate = 1
    # outputs = all
  [../]
[./Tonmergel_nomech]
    type = RedbackMaterial
    block = '2'
    gr = 1e-4
    ar = 12
    temperature = temp
    pore_pres = pore_pressure
    phi0 = 0.1
    Kc = 10
    ref_lewis_nb = 1
    ar_F = 16
    disp_x = disp_x
    ar_R = 10
    disp_y = disp_y
    eta1 = 1e3
    mu = 1
    Aphi = 1
    is_mechanics_on = true
    # total_porosity = total_porosity
    solid_thermal_expansion = 1e-3
    #gravity = '0 -9.81 0'
    solid_density = 2.3
    solid_compressibility = 1.1111111111
  [../]
  [./Tonmergel_mech]
    type = RedbackMechMaterialJ2
    block = '2'
    youngs_modulus = 10
    exponent = 5
    disp_y = disp_y
    disp_x = disp_x
    yield_stress = '0. 0.1 1. 0.1'
    temperature = temp
    pore_pres = pore_pressure
    poisson_ratio = 0.25
    ref_pe_rate = 1
    # outputs = all
 [../]
[./Malm_nomech]
    type = RedbackMaterial
    block = '3'
    gr = 1e-5
    ar = 12
    temperature = temp
    pore_pres = pore_pressure
    phi0 = 0.1
    Kc = 10
    ref_lewis_nb = 1
    ar_F = 16
    disp_x = disp_x
    ar_R = 10
    disp_y = disp_y
    eta1 = 1e3
    mu = 1
    Aphi = 1
    is_mechanics_on = true
    # total_porosity = total_porosity
    solid_thermal_expansion = 1e-3
    #gravity = '0 -9.81 0'
    solid_density = 2.68
    solid_compressibility = 1.1111111111
  [../]
  [./Malm_mech]
    type = RedbackMechMaterialJ2
    block = '3'
    youngs_modulus = 40
    exponent = 5
    disp_y = disp_y
    disp_x = disp_x
    yield_stress = '0. 1 1. 1'
    temperature = temp
    pore_pres = pore_pressure
    poisson_ratio = 0.25
    ref_pe_rate = 1
    # outputs = all
    [../]
[./Basement_nomech]
    type = RedbackMaterial
    block = '4'
    gr = 1e-5
    ar = 12
    temperature = temp
    pore_pres = pore_pressure
    phi0 = 0.1
    Kc = 10
    ref_lewis_nb = 1
    ar_F = 16
    disp_x = disp_x
    ar_R = 10
    disp_y = disp_y
    eta1 = 1e3
    mu = 1
    Aphi = 1
    is_mechanics_on = true
    # total_porosity = total_porosity
    solid_thermal_expansion = 1e-3
    #gravity = '0 -9.81 0'
    solid_density = 2.78
    solid_compressibility = 1.1111111111
  [../]
  [./Basement_mech]
    type = RedbackMechMaterialJ2
    block = '4'
    youngs_modulus = 60
    exponent = 5
    disp_y = disp_y
    disp_x = disp_x
    yield_stress = '0. 1 1. 1'
    temperature = temp
    pore_pres = pore_pressure
    poisson_ratio = 0.2
    ref_pe_rate = 1
    # outputs = all
  [../]
[]

[Constraints]
  [./slider_top]
    type = EqualValueBoundaryConstraint
    variable = disp_y
    master = 289
    slave = 4
    penalty = 1e2
  [../]
[]

[Functions]
active = 'init_gradient_T time_linear_increase pres_top'
  [./time_linear_increase]
    type = ParsedFunction
    value = 'if(t<5, 0, -1*t)'
    #vals = scale
    #vars = scale
  [../]
  [./init_gradient_T]
    type = ParsedFunction
    value = 'T_max + (y-y_min)*(T_min-T_max)/(y_max-y_min)'
    vals = '-6000            0          0          30000          28           208           0.2'
    vars = 'y_min y_max x_min x_max  T_min T_max amplitude'
  [../]
  [./init_gradient_P]
    type = ParsedFunction
    value = 'P_max + (y-y_min)*(P_min-P_max)/(y_max-y_min)'
    vals = '-6000            0          0          30000          0.01383           0.14715           0.2'
    vars = 'y_min y_max x_min x_max  P_min P_max amplitude'
  [../]
  [./pres_top]
    type = ParsedFunction
    value = 'if(t<5, 0.1*t*0.01385, 0.01385)'
  [../]
[]

[BCs]
active = 'hold_left vel_right temperature_bot temperature_top Pressure_top Pressure_bot hold_bot'
  [./hold_left]
    type = DirichletBC
    variable = disp_x
    boundary = 1
    value = 0
  [../]
  [./vel_right]
    type = FunctionDirichletBC
    variable = disp_x
    boundary = 2
    function = time_linear_increase
  [../]
  [./Pressure_top]
    type = FunctionDirichletBC
    variable = pore_pressure
    boundary = 4
    function = init_gradient_P
  [../]
  [./Pressure_bot]
    type = FunctionDirichletBC
    variable = pore_pressure
    boundary = 3
    function = init_gradient_P
  [../]
  [./hold_bot]
    type = PresetBC
    variable = disp_y
    boundary = 3
    value = 0
  [../]
  [./temperature_bot]
    type = FunctionDirichletBC
    variable = temp
    boundary = 3
    function = init_gradient_T
  [../]
  [./temperature_top]
    type = FunctionDirichletBC
    variable = temp
    boundary = 4
    function = init_gradient_T
  [../]
[]

[AuxVariables]
active = 'returnmap_iter mises_stress stress_x stress_y peeq youngs_modulus'
  [./returnmap_iter]
    order = CONSTANT
    family = MONOMIAL
    block = '1 2 3 4 5'
  [../]
  [./youngs_modulus]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./mises_stress]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./stress_x]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./stress_y]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./peeq]
    order = CONSTANT
    family = MONOMIAL
  [../]
[]

[Kernels]
active = 'temp_diff temp_td press_td press_diff'
  [./temp_td]
    type = TimeDerivative
    variable = temp
  [../]
  [./temp_diff]
    type = RedbackThermalDiffusion
    variable = temp
  [../]
  [./temp_conv]
    type = RedbackThermalConvection
    variable = temp
  [../]
  [./press_td]
    type = TimeDerivative
    variable = pore_pressure
  [../]
  [./press_diff]
    type = RedbackMassDiffusion
    variable = pore_pressure
  [../]
  [./pres_conv]
    type = RedbackMassConvection
    variable = pore_pressure
    temperature = temp
  [../]
  [./press_thermPress]
    type = RedbackThermalPressurization
    variable = pore_pressure
    temperature = temp
  [../]
  [./mh_temp]
    type = RedbackMechDissip
    variable = temp
    block= '1 2 3 4 5'
  [../]
[]

[AuxKernels]
active = 'returnmap_iter mises_stress stress_x stress_y eqv_plastic_strain'
  [./returnmap_iter]
    type = MaterialRealAux
    variable = returnmap_iter
    property = returnmap_iter
    block = '1 2 3 4 5'
  [../]
  [./mises_stress]
    type = MaterialRealAux
    variable = mises_stress
    property = mises_stress
  [../]
  [./stress_x]
    type = RankTwoAux
    rank_two_tensor = stress
    variable = stress_x
    index_i = 0
    index_j = 0
  [../]
  [./stress_y]
    type = RankTwoAux
    rank_two_tensor = stress
    variable = stress_y
    index_i = 1
    index_j = 1
  [../]
  [./eqv_plastic_strain]
    type = MaterialRealAux
    variable = peeq
    property = eqv_plastic_strain
  [../]
[]


[Postprocessors]
active = 'max_returnmap_iter'
  [./max_returnmap_iter]
    type = ElementExtremeValue
    variable = returnmap_iter
  [../]
[]

[Preconditioning]
  # active = ''
  [./SMP]
    type = SMP
    full = true
  [../]
[]

[Executioner]
  # Preconditioned JFNK (default)
  start_time = 0.0
  end_time = 6000
  dtmax = 1
  dtmin = 1e-7
  type = Transient
  num_steps = 500000
  l_max_its = 500
  solve_type = PJFNK
  petsc_options_iname = '-pc_type -pc_hypre_type -snes_linesearch_type -ksp_gmres_restart'
  petsc_options_value = 'hypre boomeramg cp 201'
  nl_abs_tol = 1e-10 # 1e-10 to begin with
  reset_dt = true
  line_search = basic
  [./TimeStepper]
    type = ReturnMapIterDT
    dt = 1e-5
    min_iter = 10
    ratio = 0.5
    max_iter = 20
    dt_max = 1e-1
    postprocessor = max_returnmap_iter
    dt_min = 1e-5
  [../]
[]

[Outputs]
  file_base = test
  output_initial = true
  exodus = true
  csv = true
  [./console]
    type = Console
    perf_log = true
    output_linear = true
    file_base = simulation
  [../]
[]

[RedbackMechAction]
  [./solid]
    disp_y = disp_y
    disp_x = disp_x
    temp = temp
  [../]
[]

[ICs]
  active = 'IC_temp IC_pressure IC_Youngs_Quartear IC_Youngs_TM IC_Youngs_Malm IC_Youngs_Base IC_Youngs_Fault'
  [./IC_temp]
    function = init_gradient_T
    variable = temp
    type = FunctionIC
  [../]
  [./IC_pressure]
    function = init_gradient_P
    variable = pore_pressure
    type = FunctionIC
  [../]
  [./IC_Youngs_Quartear]
    block = '1'
    variable = youngs_modulus
    type = FunctionNormalDistributionIC
    standard_deviation = 1.5
    maximum = 20
    minimum = 10
    mean = 15
  [../]
  [./IC_Youngs_TM]
    block = '2'
    variable = youngs_modulus
    type = FunctionNormalDistributionIC
    standard_deviation = 1
    maximum = 15
    minimum = 5
    mean = 10
  [../]
  [./IC_Youngs_Malm]
    block = '3'
    variable = youngs_modulus
    type = FunctionNormalDistributionIC
    standard_deviation = 4
    maximum = 50
    minimum = 30
    mean = 40
  [../]
  [./IC_Youngs_Base]
    block = '4'
    variable = youngs_modulus
    type = FunctionNormalDistributionIC
    standard_deviation = 6
    maximum = 70
    minimum = 50
    mean = 60
  [../]
  [./IC_Youngs_Fault]
    block = '5'
    variable = youngs_modulus
    type = FunctionNormalDistributionIC
    standard_deviation = 0.5
    maximum = 7.5
    minimum = 2.5
    mean = 5
  [../]
[]
