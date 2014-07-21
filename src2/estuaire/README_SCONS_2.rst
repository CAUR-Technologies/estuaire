Since SCons 2.xxx decided to replace the default pickle module by the cPickle module completely inside python default module dictionary, every modules that depends on feature present in pickle module but not in cPickle module (e.g. scipy) will failed miserably with different type of missing attribute error.


To ensure the original pickle modules is kept intact, you should set the environment variable SCONS_HORRIBLE_REGRESSION_TEST_HACK before calling the SConstruct.
