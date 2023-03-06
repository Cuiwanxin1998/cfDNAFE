import os
from hashlib import md5
from .utils import commonError, flatten, rmEndString
import time
import sys
import subprocess
from multiprocessing import Pool
import pandas as pd

__metaclass__ = type


class StepBase:
    def __init__(self):
        self.inputs = {}
        self.outputs = {}
        self.logpath = {}
        self.params = {}
        self.attentionSteps = ["readCounter"]

    def setInput(self, inputName, inputValue):
        if isinstance(inputName, list):
            if len(inputName) != len(inputValue):
                raise commonError("Number of input name and value not equal.")
            values = self.absolutePath(inputValue)
            for name, value in zip(inputName, values):
                self.inputs[name] = value
        else:
            self.inputs[inputName] = [self.absolutePath(inputValue)]

    # get input value by dict name (key)
    def getInput(self, inputName):
        return self.inputs[inputName]

    # get all input keys

    def getInputs(self, ):
        return list(self.inputs.keys())

    def setOutput(self, outputName, outputValue):
        if isinstance(outputName, list):
            if len(outputName) != len(outputValue):
                raise commonError("Number of output key name and value not equal.")
            values = self.absolutePath(outputValue)
            for name, value in zip(outputName, values):
                self.outputs[name] = value
        else:
            self.outputs[outputName] = self.absolutePath(outputValue)

        # get output value by dict name (key)

    def getOutput(self, outputName):
        return self.outputs[outputName]

        # get all output keys

    def getOutputs(self, ):
        return list(self.outputs.keys())


    def absolutePath(self, pathOrPathList):
        if pathOrPathList is None:
            return None
        elif isinstance(pathOrPathList, list):
            return [os.path.abspath(s) for s in pathOrPathList]
        else:
            return os.path.abspath(pathOrPathList)

    def setParam(self, paramName, paramValue):
        self.params[paramName] = paramValue

        # get input and output parameters

    def getParam(self, paramName):
        return self.params[paramName]

        # get all parameter keys

    def getParams(self, ):
        return list(self.params.keys())


    def cmdCreate(self, cmdlist):
        if not isinstance(cmdlist, list):
            raise commonError("Parameter 'cmdlist' must be a list!")
        else:
            tmp_cmd = []
            for tmp_value in cmdlist:
                if isinstance(tmp_value, dict):
                    for k, v in tmp_value.items():
                        if isinstance(v, bool) and v:
                            tmp_cmd.append(k)
                        elif isinstance(v, bool) and (not v):
                            pass
                        else:
                            tmp_cmd.append(k)
                            tmp_cmd.append(v)
                elif isinstance(tmp_value, list):
                    for l in tmp_value:
                        tmp_cmd.append(l)
                else:
                    tmp_cmd.append(tmp_value)

        cmd = " ".join([str(x) for x in tmp_cmd])
        return cmd

    # run the command line, this function is designed for single threads

    def run(self, cmds):
        if isinstance(cmds, list):#cmd is a list
            for idx, cmd in enumerate(cmds):
                print("Now, running command: {}".format(cmd))
                proc = subprocess.Popen(
                    cmd,
                    shell=True,
                    stdout=subprocess.PIPE,
                    stderr=subprocess.PIPE,
                    universal_newlines=True,
                )
                # get cmd info
                error = proc.communicate()
                exitCode = proc.returncode
                if exitCode == 0:  # successfully finished
                    print("Finish cmd")
                else:
                    print("failed %s %s" % (error))

        else:#cmd is a string
            print("Now, running command: {}".format(cmds))
            proc = subprocess.Popen(
                cmds,
                shell=True,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                universal_newlines=True,
            )
            # get cmd info
            output, error = proc.communicate()
            exitCode = proc.returncode
            if exitCode == 0:  # successfully finished
                print("Finish cmd")
            else:
                print("failed %s %s" % (output, error))


    def funRun(self, args):
        try:
            fun = args[0]
            param = args[1]
            results = fun(*param)
            flag = True
        except Exception as e:
            results = e
            flag = False

        return results, flag

    def cmdRun(self, cmd):
        proc = subprocess.run(
            cmd,
            shell=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE
        )
        # get cmd info
        mess = proc.stderr + proc.stdout
        if proc.returncode == 0:
            return mess, True
        else:
            mess1 = (
                "\n\n\nAn Error Occured During The Following Command Line Executing.\n"
            )
            mess2 = (
                "\n         Please Stop The Program To Check The Error.         \n\n\n"
            )
            mess_out = """\n^^^{}^^^\n{}^^^\n{}^^^\n{}^^^""".format(mess1, cmd, mess, mess2)
            return mess_out, False

        # time counter for multiRun

    def track_job(self, job, time_start, update_interval=3, print_interval=300):
        """
        from stackoverflow
        job: multiprocessing job
        """
        minute_count = 0
        while job._number_left > 0:
            delta_minutes = (time.time() - time_start) // print_interval
            if delta_minutes and (delta_minutes != minute_count):
                print(
                    "{0} minutes has passed since this step start.".format(
                        delta_minutes * 5
                    )
                )
                print("Tasks remaining = {0}".format(job._number_left * job._chunksize))
                minute_count = delta_minutes

            time.sleep(update_interval)

    def multiRun(self, args, func=None, nCore=1):
        """
        This function is designed for multiCore.
        multiRun(func, args, type, nCore=1)
        {P}arameters:
            args: Parameters for multirun, [[1, 2, 3], [2, 3, 4]] or ["cmd1", "cmd2", "cmd3"].
            func: function name for multicore, None mean cmd mode.
            nCore: int, how many cores to use.
        """
        # time start
        time_start = time.time()
        p = Pool(nCore)
        print("Start multicore running, master process number: {}".format(nCore))
        print(
            "Note: some command line verbose may be blocked, the program will record them in record file."
        )
        if func is None:
            results = p.map_async(self.cmdRun, args)
        else:

            reshaped_args = [[[func, x]] for x in args]
            results = p.starmap_async(self.funRun, reshaped_args)

        # print mess
        print("Subprocesses Start running......")
        print("Waiting for all subprocesses done...")

        self.track_job(
            job=results, time_start=time_start, update_interval=3, print_interval=300
        )

        p.close()
        p.join()
        print("All subprocesses done.")



    def getMaxFileNamePrefixV1(self, file):
        final_name = os.path.splitext(os.path.basename(file))[0]
        return final_name

    def getMaxFileNamePrefixV2(self, file):
        final_name = os.path.splitext(os.path.basename(file))[0]
        final_name = rmEndString(
            final_name,
            [
                ".bed.gz",
                ".bed",
                ".bam"
            ],
        )
        return final_name
