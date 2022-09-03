import os
import sys
import subprocess

REDUCE = r"bootstrapreduce"
ARG1 = '~/f4/linear-algebra-gb/f4/test/{0}.red'
ARG2 = '-w ~/f4/linear-algebra-gb/f4/test/output.rlg'

KEYWORDS = ['Wrong Answer', '*****', 'Error']

class Report:
    def __init__(self, msg, ok):
        self.msg = msg
        self.ok = ok
    def __str__(self):
        if self.ok:
            return "Tests Passed!"
        else:
            return f"Tests Failed:\n{self.msg}"

def slice_message(msg, i):
    return msg[max(0, i - 40):min(len(msg), i + 60)]

def compare_output(test, true):

    for kw in KEYWORDS:
        if kw in test:
            i = test.index(kw)
            one = slice_message(test, i)
            msg = f"----------- {kw} -----------\n{one}\n----------------------------"
            return Report(msg, False)

    ltest = list(test)
    ltrue = list(true)

    nonsign = ['\r', '\n']

    i, j = 0, 0
    while i < len(ltest) and j < len(ltrue):
        while i < len(ltest) and ltest[i] in nonsign:
            i += 1
        while j < len(ltrue) and ltrue[j] in nonsign:
            j += 1

        if i < len(ltest) and j < len(ltrue) and ltest[i] != ltrue[j]:
            one = slice_message(test, i)
            two = slice_message(true, j)
            msg = f"----------- Expected -----------\n{two}\n----------- Found -----------\n{one}\n----------------------------"
            return Report(msg, False)
        i += 1
        j += 1

    return Report("", True)

def runtests_f4(upd):
    p = subprocess.run([REDUCE, ARG2, ARG1.format("main")],
                        capture_output=True,
                        shell=False)

    output = p.stdout.decode()
    with open("f4last.rlg", 'w') as lastfile:
        lastfile.write(output)

    if upd:
        with open("f4correct.rlg", 'w') as correctfile:
            correctfile.write(output)

    with open("f4correct.rlg", 'r') as correctfile:
        correct = correctfile.read()

    report = compare_output(output, correct)

    return report

def main(argv):

    
    upd = "-upd" in argv

    if '-h' in argv or '--help' in argv:
        print("runtests.py usage:\n\t-upd to update results")

    
    print("Running F4 tests..")
    report = runtests_f4(upd)
    print(report)
    print("------------------------")

if __name__ == '__main__':
    main(sys.argv)
