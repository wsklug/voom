function(check_typedef_uint)
	MESSAGE(STATUS "Trying to compile CheckUnit.cc")
	try_compile(UINT_TYPEDEF_NEEDED ${CMAKE_BINARY_DIR} 
		${CMAKE_BINARY_DIR}/CMakeFiles/CheckUint.cc OUTPUT_VARIABLE UINT_OUTPUT)
	#MESSAGE(STATUS ${UINT_OUTPUT})
	set(UINT_TYPEDEF_NEEDED ${UINT_TYPEDEF_NEEDED} PARENT_SCOPE)
endfunction(check_typedef_uint)