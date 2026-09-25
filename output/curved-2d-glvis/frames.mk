build/release/advection_frames: tmp/curved-visualizations/advection_frames.cpp $(RELEASE_FILTER_LIBRARY)
	$(MFEM_CXX) $(CPPFLAGS) $(CXXFLAGS_COMMON) $(RELEASE_FLAGS) $< -o $@ $(RELEASE_FILTER_LIBRARY) $(MFEM_LIBS)
