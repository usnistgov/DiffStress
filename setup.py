from distutils.core import setup
setup(name='ResidualStress',
      version='0.1',
      description='A python pack to perform the residual stress analysis',
      author='Youngung Jeong',
      author_email='youngung.jeong@gmail.com',
      packages=['RS','RS.montecarlo'],
      package_dir={
        'RS':'src',
        'RS.motecarlo':'src/montecarlo'
        }
      )
