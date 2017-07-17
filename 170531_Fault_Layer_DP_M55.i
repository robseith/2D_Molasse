# Gravity step without displacement (Stressing step)
# Mechanical layering with fault
# Materials redback material mech (J2) and nonmech
# Disp. BC: left Disp_x = 0; right Disp_x = function ; bottom Disp_y = 0
# Constrain - Slider: Top
# Temp. BC: Top = function depth / constant, Bot = function depth / constant
# IC: Temp- function depth;
# Executioner with "Returnmap" - timestepper

[GlobalParams]
  time_factor = 1.e-3
[]

[Mesh]
  type = FileMesh
  file = Layered_Fault_TVD_subdiv_Real_fine_quad.msh
  dim = 2
[]

[Variables]
  [./disp_x]
    order = FIRST
    family = LAGRANGE
  [../]
  [./disp_y]
    order = FIRST
    family = LAGRANGE
  [../]
  #[./disp_z]
  #  order = FIRST
  #  family = LAGRANGE
  #[../]
  [./temp]
  [../]
  #[./pore_pressure]
  #[../]
[]

[Materials]
  [./Fault_nomech]
    type = RedbackMaterial
    block = '5 6'
    disp_x = disp_x
    disp_y = disp_y
    #disp_z = disp_z
    #pore_pres = pore_pressure
    temperature = temp
    #total_porosity = total_porosity
    is_mechanics_on = true
    gr = 1e-7.26
    ar = 7.26
    alpha_1 = -0.256
    alpha_2 = 3.53
    delta = 1e-3
    ref_lewis_nb = 1
    pressurization_coefficient = 1e-5
    phi0 = 0.1
    Kc = 10
    eta1 = 1e3
    mu = 1
    Aphi = 1
    solid_thermal_expansion = 1e-3
    init_from_functions__function_names = gravity_ramp
    init_from_functions__params = gravity
    solid_density = 2.5
    solid_compressibility = 7.5e-8 # 1.1111111111
  [../]
  [./Fault_mech]
    type = RedbackMechMaterialDP
    block = '5 6'
    disp_x = disp_x
    disp_y = disp_y
    #disp_z = disp_z
    temperature = temp
    #pore_pres = pore_pressure
    #total_porosity = total_porosity
    yield_stress = '0. 5e4 1. 5e4'
    slope_yield_surface = 0
    youngs_modulus = 20e6
    exponent = 1
    poisson_ratio = 0.25
  [../]
  [./Quartaer_nomech]
    type = RedbackMaterial
    block = 1
    disp_x = disp_x
    disp_y = disp_y
    #disp_z = disp_z
    #pore_pres = pore_pressure
    temperature = temp
    #total_porosity = total_porosity
    is_mechanics_on = true
    gr = 1e-5
    ar = 12
    ar_F = 12
    ar_R = 12
    phi0 = 0.1
    Kc = 10
    eta1 = 1e3
    mu = 1
    Aphi = 1
    solid_thermal_expansion = 1e-3
    init_from_functions__function_names = gravity_ramp
    init_from_functions__params = gravity
    solid_density = 2.35
    solid_compressibility = 1.26e-7 # 1.1111111111
  [../]
  [./Quartaer_mech]
    type = RedbackMechMaterialElastic
    block = 1
    disp_x = disp_x
    disp_y = disp_y
    #disp_z = disp_z
    temperature = temp
    #pore_pres = pore_pressure
    #total_porosity = total_porosity
    youngs_modulus = 10e6
    exponent = 4
    poisson_ratio = 0.3
  [../]
  [./Tonmergel_nomech]
    type = RedbackMaterial
    block = 2
    disp_x = disp_x
    disp_y = disp_y
    #disp_z = disp_z
    #pore_pres = pore_pressure
    temperature = temp
    #total_porosity = total_porosity
    is_mechanics_on = true
    gr = 1e-7.26
    ar = 7.26
    alpha_1 = -0.256
    alpha_2 = 3.53
    delta = 1e-3
    ref_lewis_nb = 1
    pressurization_coefficient = 1e-5
    phi0 = 0.1
    Kc = 10
    eta1 = 1e3
    mu = 1
    Aphi = 1
    solid_thermal_expansion = 1e-3
    init_from_functions__function_names = gravity_ramp
    init_from_functions__params = gravity
    solid_density = 2.3
    solid_compressibility = 7.5e-8 # 1.1111111111
  [../]
  [./Tonmergel_mech]
    type = RedbackMechMaterialDP
    block = 2
    disp_x = disp_x
    disp_y = disp_y
    #disp_z = disp_z
    temperature = temp
    #pore_pres = pore_pressure
    #total_porosity = total_porosity
    yield_stress = '0. 4.5e4 1. 4.5e4'
    slope_yield_surface = 0
    youngs_modulus = 20e6
    exponent = 5
    poisson_ratio = 0.25
  [../]
  [./Malm_nomech]
    type = RedbackMaterial
    block = '3'
    disp_x = disp_x
    disp_y = disp_y
    #disp_z = disp_z
    temperature = temp
    #pore_pres = pore_pressure
    #total_porosity = total_porosity
    is_mechanics_on = true
    gr = 1e-8
    ar = 8
    alpha_1 = 4.95
    alpha_2 = 0.85
    delta = 1e-3
    ref_lewis_nb = 1
    pressurization_coefficient = 1e-5
    phi0 = 0.1
    Kc = 10
    eta1 = 1e3
    mu = 1
    Aphi = 1
    solid_thermal_expansion = 1e-3
    init_from_functions__function_names = gravity_ramp
    init_from_functions__params = gravity
    solid_density = 2.65
    solid_compressibility = 5e-8 # 1.1111111111
  [../]
  [./Malm_mech]
    # outputs = all
    type = RedbackMechMaterialDP
    block = '3'
    disp_x = disp_x
    disp_y = disp_y
    #disp_z = disp_z
    temperature = temp
    #pore_pres = pore_pressure
    #total_porosity = total_porosity
    yield_stress = '0. 5.5e4 1. 5.5e4'
    slope_yield_surface = 0
    youngs_modulus = 30e6
    exponent = 2
    poisson_ratio = 0.25
  [../]
  [./Base_nomech]
    type = RedbackMaterial
    block = '4'
    disp_x = disp_x
    disp_y = disp_y
    #disp_z = disp_z
    #pore_pres = pore_pressure
    temperature = temp
    #total_porosity = total_porosity
    is_mechanics_on = true
    gr = 1e-8
    ar = 8
    alpha_1 = 4.95
    alpha_2 = 0.85
    delta = 1e-3
    ref_lewis_nb = 1
    pressurization_coefficient = 1e-5
    phi0 = 0.1
    Kc = 10
    eta1 = 1e3
    mu = 1
    Aphi = 1
    solid_thermal_expansion = 1e-3
    init_from_functions__function_names = gravity_ramp
    init_from_functions__params = gravity
    solid_density = 2.780
    solid_compressibility = 3.6e-8 # 1.1111111111
  [../]
  [./Base_mech]
    # outputs = all
    type = RedbackMechMaterialDP
    block = '4'
    disp_x = disp_x
    disp_y = disp_y
    #disp_z = disp_z
    temperature = temp
    #pore_pres = pore_pressure
    #total_porosity = total_porosity
    slope_yield_surface = 0
    yield_stress = '0. 1.5e5 1. 1.5e5'
    youngs_modulus = 50e6
    exponent = 1
    poisson_ratio = 0.2
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
  [./init_gradient_T]
    type = ParsedFunction
    value = 'T_max + (y-y_min)*(T_min-T_max)/(y_max-y_min)'
    vals = '-6500            0          0          30000          10           208           0.2'
    vars = 'y_min y_max x_min x_max  T_min T_max amplitude'
  [../]
  [./gravity_ramp]
    type = ParsedVectorFunction
    value_x = 0
    value_y = -9.81*tanh(0.001*t*t)
    value_z = 0
  [../]
  [./time_linear_increase]
    type = ParsedFunction
    value = 'if(t<20, 0, 0.001*tanh((t-20)*0.1)*(t-20))' # Deformation rate of 1 mm/a ~ 0.001 m/a
    #vals = scale
    #vars = scale
  [../]
[]

[BCs]
  active = 'vel_left hold_bot hold_right temperature_bot temperature_top'
  [./hold_left]
    type = PresetBC
    variable = disp_x
    boundary = 1
    value = 0
  [../]
  [./vel_left]
    type = FunctionDirichletBC
    variable = disp_x
    boundary = 1
    function = time_linear_increase
  [../]
  [./hold_right]
    type = PresetBC
    variable = disp_x
    boundary = 2
    value = 0
  [../]
  [./vel_right]
    type = FunctionDirichletBC
    variable = disp_x
    boundary = 2
    function = time_linear_increase
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
  [./returnmap_iter]
    order = CONSTANT
    family = MONOMIAL
    block = '1 2 3 4 5 6'
  [../]
  [./youngs_modulus]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./mises_stress]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./mech_diss]
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
  [./stress_z]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./stress_xy]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./stress_yz]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./stress_xz]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./peeq]
    order = CONSTANT
    family = MONOMIAL
  [../]
  #[./stress]
  #  order = CONSTANT
  #  family = MONOMIAL
  #[../]
  #[./total_porosity]
  #  order = FIRST
  #  family = MONOMIAL
  #[../]
  #[./mech_porosity]
  #  order = FIRST
  #  family = MONOMIAL
  #[../]
[]

[Kernels]
  active = 'temp_diff temp_dissip td_temp temp_conv'#press_diff PoroMech press_thermPress td_press
  [./td_press]
    type = TimeDerivative
    variable = pore_pressure
  [../]
  [./PoroMech]
    type = RedbackPoromechanics
    variable = pore_pressure
  [../]
  [./press_thermPress]
    type = RedbackThermalPressurization
    variable = pore_pressure
    temperature = temp
  [../]
  [./temp_diff]
    type = RedbackThermalDiffusion
    variable = temp
  [../]
  [./temp_dissip]
   type = RedbackMechDissip
   variable = temp
  [../]
  [./td_temp]
    type = TimeDerivative
    variable = temp
  [../]
  [./press_diff]
    type = RedbackMassDiffusion
    variable = pore_pressure
  [../]
  [./temp_conv]
    type = RedbackThermalConvection
    variable = temp
  [../]
  [./pres_conv]
    type = RedbackMassConvection
    variable = pore_pressure
    temperature = temp
  [../]
[]

[AuxKernels]
  [./returnmap_iter]
    type = MaterialRealAux
    variable = returnmap_iter
    property = returnmap_iter
    block = '1 2 3 4 5 6'
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
  [./stress_z]
    type = RankTwoAux
    rank_two_tensor = stress
    variable = stress_z
    index_i = 2
    index_j = 2
  [../]
  [./stress_xy]
    type = RankTwoAux
    rank_two_tensor = stress
    variable = stress_xy
    index_i = 0
    index_j = 1
  [../]
  [./stress_yz]
    type = RankTwoAux
    rank_two_tensor = stress
    variable = stress_yz
    index_i = 1
    index_j = 2
  [../]
  [./stress_xz]
    type = RankTwoAux
    rank_two_tensor = stress
    variable = stress_xz
    index_i = 0
    index_j = 2
  [../]
  #[./stress]
  #  type =  MaterialRealTensorValueAux
  #  variable = stress
  #  property = rank_two_tensor
    #rank_two_tensor = stress
    #variable = L2norm
    #scalar_type = L2norm
  #[../]
  [./eqv_plastic_strain]
    type = MaterialRealAux
    variable = peeq
    property = eqv_plastic_strain
  [../]
  [./mech_dissipation]
    type = MaterialRealAux
    variable = mech_diss
    property = mechanical_dissipation_mech
  [../]
  #[./total_porosity]
  #  type = RedbackTotalPorosityAux
  #  variable = total_porosity
  #  mechanical_porosity = mech_porosity
  #[../]
  #[./mech_porosity]
  #  type = MaterialRealAux
  #  variable = mech_porosity
  #  execute_on = timestep_end
  #  property = mechanical_porosity
  #[../]
[]

[Postprocessors]
  [./evaluations]           # Cumulative residual calculations for simulation
    type = NumResidualEvaluations
  [../]
  [./active_time]           # Time computer spent on simulation
    type = RunTime
    time_type = active
  [../]
  [./max_returnmap_iter]
    type = ElementExtremeValue
    variable = returnmap_iter
  [../]
  [./dt]
    type = TimestepSize
  [../]
  [./TM_MisS_mid]
    type = PointValue
    point = '13000 -3000 0'
    variable = mises_stress
  [../]
  [./TM_MisS_base]
    type = PointValue
    point = '13000 -3349 0'
    variable = mises_stress
  [../]
  [./TM_SX_mid]
    type = PointValue
    point = '13000 -3000 0'
    variable = stress_x
  [../]
  [./TM_SY_mid]
    type = PointValue
    point = '13000 -3000 0'
    variable = stress_y
  [../]
  [./Malm_MisS_mid]
    type = PointValue
    point = '13000 -3750 0'
    variable = mises_stress
  [../]
  [./Malm_MisS_base]
    type = PointValue
    point = '13000 -4146 0'
    variable = mises_stress
  [../]
  [./Malm_SX_mid]
    type = PointValue
    point = '13000 -3750 0'
    variable = stress_x
  [../]
  [./Malm_SY_mid]
    type = PointValue
    point = '13000 -3750 0'
    variable = stress_y
  [../]
  [./Base_MisS_mid]
    type = PointValue
    point = '13000 -5000 0'
    variable = mises_stress
  [../]
  [./Base_MisS_base]
    type = PointValue
    point = '13000 -6500 0'
    variable = mises_stress
  [../]
  [./Base_SX_mid]
    type = PointValue
    point = '13000 -5000 0'
    variable = stress_x
  [../]
  [./Base_SY_mid]
    type = PointValue
    point = '13000 -5000 0'
    variable = stress_y
  [../]
  [./Disp_x_left]
    type = PointValue
    point = '0 -5000 0'
    variable = disp_x
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
  end_time = 2e5 # ~ 1.6% strain (def endtime 1 Mio years with deform. rate of 1 mm/a ~3.3% Strain)
  type = Transient
  num_steps = 500000
  l_max_its = 500
  solve_type = PJFNK
  petsc_options_iname = '-pc_type -pc_hypre_type -snes_linesearch_type -ksp_gmres_restart'
  petsc_options_value = 'hypre boomeramg cp 201'
  nl_abs_tol = 1e-2
  reset_dt = true
  line_search = basic
  [./TimeStepper]
    type = ReturnMapIterDT
    dt = 1e-9
    min_iter = 10
    ratio = 0.7
    max_iter = 20
    dt_max = 1000
    postprocessor = max_returnmap_iter
    dt_min = 1e-12
  [../]
[]

[Debug]
  show_var_residual_norms = true
[]

[Outputs]
  file_base = 170531_extension_70_M55
  output_initial = true
  exodus = true
  csv = true
  [./console]
    type = Console
    perf_log = true
    output_linear = true
    file_base = simulation
  [../]
  [./exodus_output]
    file_base = 170531_Fault_Layer_SI_DP_70_M55
    elemental_as_nodal = true
    #output_material_properties = true
    execute_nodal_on = 'initial timestep_end'
    execute_vector_postprocessors_on = 'initial timestep_end'
    type = Exodus
  [../]
  #[./CSV_output]
  #  file_base = 170331_Fault_Layer_SI_DP_70_M46
  #  elemental_as_nodal = true
  #  #output_material_properties = true
  #  execute_nodal_on = 'final' #'initial timestep_end'
  #  type = CSV
  #[../]
[]
#[Outputs]
#  [./out]
#    file_base = 170331_test
#    type = Exodus
#    output_material_properties = true
#    show_material_properties = 'tensor_properties'
#  [../]
#[]

[RedbackMechAction]
  [./solid]
    disp_y = disp_y
    disp_x = disp_x
    #disp_z = disp_z
    temp = temp
    #pore_pres = pore_pressure
  [../]
[]

[ICs]
  #  active = 'IC_temp'
  [./IC_temp]
    function = init_gradient_T
    variable = temp
    type = FunctionIC
  [../]
  #[./press_ic]
  #  variable = pore_pressure
  #  type = ConstantIC
  #  value = 0
  #[../]
  [./IC_Youngs_Quartear]
    block = '1'
    variable = youngs_modulus
    type = FunctionNormalDistributionIC
    standard_deviation = 1e6
    maximum = 15e6
    minimum = 5e6
    mean = 10e6
  [../]
  [./IC_Youngs_TM]
    block = '2'
    variable = youngs_modulus
    type = FunctionNormalDistributionIC
    standard_deviation = 2e6
    maximum = 25e6
    minimum = 15e6
    mean = 20e6
  [../]
  [./IC_Youngs_Malm]
    block = '3'
    variable = youngs_modulus
    type = FunctionNormalDistributionIC
    standard_deviation = 3e6
    maximum = 35e6
    minimum = 25e6
    mean = 30e6
  [../]
  [./IC_Youngs_Base]
    block = '4'
    variable = youngs_modulus
    type = FunctionNormalDistributionIC
    standard_deviation = 5e6
    maximum = 55e6
    minimum = 45e6
    mean = 50e6
  [../]
  [./IC_Youngs_Fault_Malm]
    block = '5 6'
    variable = youngs_modulus
    type = FunctionNormalDistributionIC
    standard_deviation = 2e6
    maximum = 25e6
    minimum = 15e6
    mean = 20e6
  [../]
[]
