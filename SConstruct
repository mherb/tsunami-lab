import os

env = Environment(tools=['default','cxxtest'])

# source and test file paths
src_files = Glob('src/*.cpp')
tests_files = Glob('tests/*.h')

# build target dir
buildDir = '#build'

env.Append(CPPPATH=['#src'])

# Unit tests
env['CXXTEST_SKIP_ERRORS'] = True
tests = env.CxxTest(os.path.join(buildDir, 'testrunner'), tests_files + src_files)

# disable default target
env.Default('check')
