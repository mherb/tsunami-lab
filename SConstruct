import os

env = Environment(tools=['default','cxxtest'])

srcDir = '#src'
testsDir = '#tests'
buildDir = '#build'

srcFiles = Glob(os.path.join(srcDir, '*.cpp'))
testsFiles = Glob(os.path.join(testsDir, '*.h'))

env.Append(CPPPATH=[srcDir])

# Unit tests
env['CXXTEST_SKIP_ERRORS'] = True
tests = env.CxxTest(os.path.join(buildDir, 'testrunner'), testsFiles + srcFiles)

env.Default('check')
