module cam_debug

  implicit none
  private
  save

  logical, parameter, public :: l_debug = .False. 
  logical, parameter, public :: l_summary_debug = .False. 
  logical, parameter, public :: l_masscon_debug = .False. 
  logical, parameter, public :: l_step2_debug = .False. 

end module cam_debug
