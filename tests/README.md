# Intensiphy Testing Suite

Welcome to Intensiphy's Testing Suite. The tests here will help you check Intensiphy's functionality if you suspect the program isn't working correctly. You can run each test individually or you can run them all together (recommended). You can run these tests as many times as you like as the output directories produced by each test run will be removed prior to new tests being run.

## Quick Start

At the command line, change directories into the `/tests` directory.

```
cd [PATH/TO]/Intensiphy/tests
```

Run this command to test each major aspect of Extensiphy at once.

```
pytest ./unit_tests.py
```

Thats it. Once the tests are complete, you should hopefully see this output:

```
unit_tests.py ...                                                                                                                        [100%]

========================================================= 3 passed, 1 warning in 0.31s =========================================================

```

If any of these tests have failed and show a `FAILED` result, there is a problem with Intensiphy  
To get help with your installation of Intensiphy, contact: jtoscanifield@ucmerced.edu
