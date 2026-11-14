r""" Shared harness for the repo's script-style test files.

    from gridsearch.testing import check, report
    check('name of the assertion', condition, 'optional extra detail')
    report()   # prints the verdict and exits nonzero on any failure

Importing this module reconfigures stdout/stderr to UTF-8. That is not cosmetic: every test file prints
Greek parameter names, and on Windows a redirected stdout defaults to the ANSI codepage, so the first `β`
raises UnicodeEncodeError and a passing suite reports FAIL for a purely clerical reason. Do not remove it,
and do not add a module that prints Greek without importing from here.

State is module-level, so a test file that imports another test file's module shares one tally -- which is
what the model suites want (`import test as testmod` pulls in the shared instance, not assertions).
"""
import sys


def utf8Stdout():
    """ Make stdout/stderr UTF-8 whatever they are attached to. Called on import; also called explicitly
    by `python/runTests.py`, which prints Greek suite names before any suite runs. """
    for stream in (sys.stdout, sys.stderr):
        try:
            stream.reconfigure(encoding = 'utf-8')
        except (AttributeError, ValueError):    # not a TextIOWrapper (captured/piped by a harness)
            pass


utf8Stdout()

_passed = 0
_failed = 0


def check(name, cond, extra = ''):
    """ Record one assertion. `cond` is coerced with bool(), so numpy scalars are fine but arrays are not
    -- wrap those in .all()/np.allclose at the call site, where the intent is visible. """
    global _passed, _failed
    cond = bool(cond)
    if cond:
        _passed += 1
    else:
        _failed += 1
    print(('PASS' if cond else 'FAIL') + '  ' + name + ('  ' + extra if extra else ''))
    return cond


def passed():
    """ True while no check has failed. """
    return _failed == 0


def report(exit = True):
    """ Print the verdict; exit nonzero on failure unless `exit = False`. The counts are what the runner
    (`python/runTests.py`) parses, so keep the last line's shape. """
    print()
    print('{}  ({} passed, {} failed)'.format('ALL PASS' if _failed == 0 else 'FAILURES ABOVE',
                                              _passed, _failed))
    if exit:
        sys.exit(0 if _failed == 0 else 1)
    return _failed == 0


def reset():
    """ Clear the tally. Only for a runner that executes several suites in one process. """
    global _passed, _failed
    _passed = _failed = 0
