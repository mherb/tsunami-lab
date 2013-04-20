import os

env = Environment(tools=['default','cxxtest'])

srcDir = '#src'
testsDir = '#tests'
buildDir = '#build'

srcFiles = Glob(os.path.join(srcDir, '*.cpp'))
testsFiles = Glob(os.path.join(testsDir, '*.h'))

env.Append(CPPPATH=[srcDir])

# Pass optimize=<level> to set optimization level
# by default level 0 (no optimzation) is used
if ARGUMENTS.get('optimize', 0):
    env.Append(CCFLAGS=['-O' + ARGUMENTS.get('optimize', '0')])

# Pass debug=1 for a debug build
if ARGUMENTS.get('debug', 0):
    env.Append(CCFLAGS=['-g'])

# Pass compiler=<name> to use <name> as compiler instead of g++
if ARGUMENTS.get('compiler', 0):
    env['CXX'] = ARGUMENTS.get('compiler')

# Unit tests
env['CXXTEST_SKIP_ERRORS'] = True
tests = env.CxxTest(os.path.join(buildDir, 'testrunner'), testsFiles + srcFiles)

env.Default(tests)
