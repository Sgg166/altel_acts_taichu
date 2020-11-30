# project(millepede LANGUAGES Fortran C)

enable_language(Fortran)
enable_language(C)

set(PEDE_SRC 
  Dbandmatrix.f90
  linesrch.f90
  mille.f90
  minresDataModule.f90
  minresModule.f90
  minresqlpBlasModule.f90
  minresqlpDataModule.f90
  minresqlpModule.f90
  mpbits.f90
  mpdalc.f90
  mpdef.f90
  mphistab.f90
  mpmod.f90
  mpnum.f90
  mpqldec.f90
  mptest1.f90
  mptest2.f90
  mptext.f90
  pede.f90
  randoms.f90
  vertpr.f90
  )

add_library(objpede OBJECT ${PEDE_SRC})
target_compile_options(objpede PRIVATE "-DREAD_C_FILES")  
set_target_properties(objpede PROPERTIES COMPILE_FLAGS "-fautomatic -fno-backslash -cpp")

add_library(objpede_fmain OBJECT pede_main.f90)
target_compile_options(objpede_fmain PRIVATE "-DREAD_C_FILES")
set_target_properties(objpede_fmain PROPERTIES COMPILE_FLAGS "-fautomatic -fno-backslash -cpp")

# find_package(OpenMP)
# if(OpenMP_Fortran_FOUND)
#   get_property(objpede_flags TARGET objpede PROPERTY COMPILE_FLAGS)
#   set(objpede_flags "${objpede_flags} ${OpenMP_Fortran_FLAGS}")
#   set_target_properties(objpede PROPERTIES COMPILE_FLAGS ${objpede_flags})
#   set_target_properties(objpede_fmain PROPERTIES COMPILE_FLAGS ${objpede_flags})
# endif()

add_library(objreadc OBJECT readc.c)
target_include_directories(objreadc PRIVATE $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>)
target_compile_options(objreadc PRIVATE "-Df2cFortran")

# add_library(pedelib SHARED $<TARGET_OBJECTS:objpede> $<TARGET_OBJECTS:objreadc>)
# if(OpenMP_Fortran_FOUND)
#   target_link_libraries(pedelib OpenMP::OpenMP_Fortran)
# endif()

add_executable(pede $<TARGET_OBJECTS:objpede_fmain> $<TARGET_OBJECTS:objpede> $<TARGET_OBJECTS:objreadc>)
add_executable(pede_cmain pede_main.cc $<TARGET_OBJECTS:objpede> $<TARGET_OBJECTS:objreadc>) #testing

install(TARGETS pede pede_cmain
  EXPORT ${PROJECT_NAME}Targets
  RUNTIME       DESTINATION bin      COMPONENT runtime
  LIBRARY       DESTINATION lib      COMPONENT runtime
  ARCHIVE       DESTINATION lib      COMPONENT devel
  PUBLIC_HEADER DESTINATION include  COMPONENT devel
  RESOURCE      DESTINATION resource COMPONENT runtime
  )
