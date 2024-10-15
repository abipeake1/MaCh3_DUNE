set(APPS EventRatesDUNE EventRates_DUNE_shiftingparameters SigmaVariationsDUNE LLHScansDUNE jointFitDUNE)

foreach(a ${APPS})
	  add_executable(${a} ${a}.cpp)
	  target_link_libraries(${a} SamplePDFDUNE)
	  install(TARGETS ${a} DESTINATION bin)
endforeach()
