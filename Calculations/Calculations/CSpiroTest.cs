using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Diagnostics;
using System.IO;
using System.Resources;

namespace Calculations
{
    public class CSpiroTest
    {
        //Variables
        public int MVersion;
        public DateTime MTimeAcquired;

        public CSpireMeas MTestBarometricPressure;
        public CSpireMeas MTestTemperature;

        public uint MSensorId;
        public bool MSensorCalibrated;
        public double MSensorCorrFactor;
        public DateTime MCalibrationTime;
        public bool MUseDryAirConditions;

        public bool MTestAccepted;
        public bool MTestAdequate;
        public string MStrReasonNotAdequate;
        public bool MIsBestTest;

        public short MPressStartVal;
        //See if I can change CByteArray into a short array for output
        //Currently making class for CByteArray to code CByteArray MPressureData;
        public  List<Byte> MPressureData = new List<Byte>();
        public short MTempStartVal;
        public short MPrevPressureAdValue;
        public short MPrevTempAdValue;

        public bool MTestStarted;
        public  int MStartTestIndex;
        public bool MSensorZeroed;
        public int MSensorZeroValue;
        public int MNumSamplesReceived;
        public int MAbsoluteStartOfExp;

        public bool MUseInspOrExpirForStart;
        public  double MVolumeSum;

        protected CIirFilter MpIirFilter;

        public  Globals.ConvTableType MConvTable = new Globals.ConvTableType(Globals.DefaultEntryLength());

        public double MCalVol;
        public bool MInitialFlowDetected;
        public int MNumPtsAboveMinFlow;
        public double MPositiveFlowVolume;

        public  int MSampleRate;

        public  double MPrevVolValue;
        public  double MPrevFlowValue;
        public  double MVolVsTimeStartVal;
        public  double MFlowVsTimeStartVal;

        public int MDemoModeSampleCount;

        public int MSensorZeroMinValue;
        public int MSensorZeroMaxValue;

        public bool MSelectedByOperator;
        public bool MInhalationStarted;

        //CWordArray = string List in C#
        public  List<ushort> MVolVsTimeWords = new List<ushort>();
        public  List<ushort> MFlowVsTimeWords = new List<ushort>();

        //Converted the next 2 variables from CArray<double, double> to a list of structs of thype double,double.
        //Using 2 different structs where one is flow and time while the other is volume and time.
        //Doing above statement for clarity, but I suggest it be changed when code is to be refactored
        //Same with in Global.cs Statement above PRESS_FLOW_TYPE
        //8/30/2017 Update: Gongyu said that CArray <double,double> is just a list of type double

        public  List<double> MFlowVsTime = new List<double>();
        public  List<double> MVolVsTime = new List<double>();

        public int MStartOfInspiration;

        public uint MFvcTestStatus;

        //Constructor
        public CSpiroTest()
        {
            MVersion = Globals.CSpiroTestVersion;
            MTimeAcquired = DateTime.Now;

            MSensorId = 0;
            MSensorCalibrated = false;
            MSensorCorrFactor = 1.0;
            MCalibrationTime = DateTime.Now;
            MUseDryAirConditions = false;

            MTestAccepted = false;
            MTestAdequate = false;
            MStrReasonNotAdequate = string.Empty;
            MIsBestTest = false;

            MPressStartVal = 0;
            MTempStartVal = 0;
            MPrevPressureAdValue = 0;
            MPrevTempAdValue = 0;

            MTestStarted = false;
            MStartTestIndex = 0;
            MSensorZeroed = false;
            MSensorZeroValue = 0;
            MNumSamplesReceived = 0;
            MAbsoluteStartOfExp = 0;

            MUseInspOrExpirForStart = false;   // Default to expiration for test start
            MVolumeSum = 0.0;

            MTestBarometricPressure = new CSpireMeas();
            MTestTemperature = new CSpireMeas();
            MTestBarometricPressure.MMeasVal = Globals.DefaultBarometicPressure;
            MTestBarometricPressure.MIsValid = true;
            MpIirFilter = new CIirFilter(Globals.SpiroFilt);

            for (var i = 0; i < Globals.DefaultEntryLength(); i++)
            {
                MConvTable.PressFlow[i].FlowVal = Globals.DefConvTable[i].FlowVal;
                MConvTable.PressFlow[i].PressVal = Globals.DefConvTable[i].PresVal;
            }

            MCalVol = 0.0;
            MInitialFlowDetected = false;
            MNumPtsAboveMinFlow = 0;
            MPositiveFlowVolume = 0.0;

            MSampleRate = Globals.SpiroSampleRate;

            MPrevVolValue = 0.0;
            MPrevFlowValue = 0.0;
            MVolVsTimeStartVal = 0.0;
            MFlowVsTimeStartVal = 0.0;

            MDemoModeSampleCount = 0;

            MSensorZeroMinValue = 65536;
            MSensorZeroMaxValue = -65536;

            MSelectedByOperator = false;
            MInhalationStarted = false;

        }

        //Functions

        //Gongyu said that this function does not need to be completed because no one references this function
        //public void ReadConvTableFromFile(String strFileName)
        //{
        //    string f;
        //    String strLine, strPressure, strFlow;
        //    int TabIndex;
        //    int i = 0;

        //    if (File.Exists(strFileName))
        //    {
        //        f = System.IO.File.ReadAllText(strFileName);

        //    }
        //    else
        //    {
        //        return;
        //    }
        //}

        //There is a serialize function that goes here
        public virtual void Serialize()
        {
            throw new NotImplementedException("Not Implemented because we need to use old c++ code for legacy deserialization");

        }


        //The reason these functions was done was because it had to be stored
        //in bytes by using the first difference method.

        public virtual void InitializeStorage(short firstAdPressureValue, short firstAdTempValue)
        {
            MPressStartVal = firstAdPressureValue;
            MTempStartVal = firstAdTempValue;
            MPrevPressureAdValue = firstAdPressureValue;
            MPrevTempAdValue = firstAdTempValue;

            MPressureData.Clear();

        }
        public virtual void StoreSamples(short adPressureValue, short adTempValue)
        {
            var tempNewPressValue = (short)(adPressureValue - MPrevPressureAdValue);
            var tempNewTempValue = (short)(adTempValue - MPrevTempAdValue);
            MPressureData.Add((Byte)tempNewPressValue);
            MPrevPressureAdValue = adPressureValue;
            MPrevTempAdValue = adTempValue;
        }

        //These functions are used for version 2 or later
        public void InitializeStorage(double firstVolVsTimeValue, double firstFlowVsTimeValue)
        {
            MPrevVolValue = firstVolVsTimeValue;
            MPrevFlowValue = firstFlowVsTimeValue;
            MVolVsTimeStartVal = firstVolVsTimeValue;
            MFlowVsTimeStartVal = firstFlowVsTimeValue;

            MVolVsTimeWords.Clear();
            MFlowVsTimeWords.Clear();
        }

        public void StoreSamples(double volVsTimeValue, double flowVsTimeValue)
        {
            short tempVolDelta;
            short tempFlowDelta;

            if (volVsTimeValue - MPrevVolValue >= 0.0)
            {
                tempVolDelta = (short)((volVsTimeValue - MPrevVolValue + 0.0005) * 1000);
            }
            else
            {
                tempVolDelta = (short)((volVsTimeValue - MPrevVolValue - 0.0005) * 1000);
            }

            if (flowVsTimeValue - MPrevFlowValue >= 0)
            {
                tempFlowDelta = (short)((flowVsTimeValue - MPrevFlowValue + 0.0005) * 1000);
            }
            else
            {
                tempFlowDelta = (short)((flowVsTimeValue - MPrevFlowValue - 0.0005) * 1000);
            }

            MVolVsTimeWords.Add((ushort)tempVolDelta);
            MFlowVsTimeWords.Add((ushort)tempFlowDelta);
            MPrevVolValue = volVsTimeValue;
            MPrevFlowValue = flowVsTimeValue;
        }
        //  Copied this, had void as a parameter
        //  Each derived class must provide a version of this function to 
        //  perform the measurements for the appropriate test.
        //  1 Jan 97 PHR

        public virtual void CalculateMeasurements()
        {

        }
        public  double AdToFlow(short adValue, short adValueZero, bool useCorrectionFactor)
        {
            var i = 0;
            var done = false;
            var upperLimit = MConvTable.NumEntries - 1;
            var samplePoint = (short)(adValue - adValueZero);
            var pressureVal = samplePoint * Globals.InchesPerBit;
            var flowVal = 0.0;
            var signCorrection = 1.0;


            if ((adValue - adValueZero) == 0)
            {
                return (0.0);
            }
            //This code never executes because MConvTable is never cleared/<= 0
            //Must be for error checking
            if (MConvTable.NumEntries <= 0)
            {
                return (0.0);
            }

            //Both if statements could never execute, unless MConvTable changes in some other function
            if (pressureVal <= MConvTable.PressFlow[0].PressVal)
            {
                return (MConvTable.PressFlow[0].FlowVal * signCorrection);
            }
            if (pressureVal > MConvTable.PressFlow[MConvTable.NumEntries - 1].PressVal)
            {
                return (MConvTable.PressFlow[MConvTable.NumEntries - 1].FlowVal * signCorrection);
            }

            while (i < upperLimit && !done)
            {
                if (pressureVal >= MConvTable.PressFlow[i].PressVal && pressureVal <= MConvTable.PressFlow[i + 1].PressVal)
                {
                    flowVal = MConvTable.PressFlow[i].FlowVal + ((pressureVal - MConvTable.PressFlow[i].PressVal) * (MConvTable.PressFlow[i + 1].FlowVal - MConvTable.PressFlow[i].FlowVal)) /
                    (MConvTable.PressFlow[i + 1].PressVal - MConvTable.PressFlow[i].PressVal);
                    done = true;
                }
                i++;
            }

            if (Math.Abs(flowVal) < 0.02)
            {
                flowVal = 0; // 0.05,07 Feb, 2006, remove small drift
            }
            return (flowVal * signCorrection);
        }

        //Equation 2.3
        public  double VMATPS_To_VMBTPS(double vmatpd, double pb, double tempC)
        {
            double vmbtps;
            vmbtps = vmatpd * ((pb - Ph20(tempC)) / (pb - 47.0) * (310.0 / (273.0 + tempC)));
            return (vmbtps);
        }

        //This function always returns a constant of 37.73498 after the functoin has been updated
        //Equation 2.4
        public  double Ph20(double tempC)
        {
            double dPh20;
            //    dPH20 = 47.07 * pow(10.0, 6.36*(TempC-37.0)/(232.0+TempC));
            // 6 Jul 00 PHR
            dPh20 = 47.07 * Math.Pow(10.0, 6.36 * (33.0 - 37.0) / (232.0 + 33.0));
            return (dPh20);
        }

        //Equation is similar to 2.3, where the different is pb is either minused by ph20 or isn't
        //Equation 2.5
        public  double VMATPD_To_VMBTPS(double vmatpd, double pb, double tempC)
        {
            double vmbtps;
            vmbtps = vmatpd * (pb / (pb - 47.0) * 310.0 / (273.0 + tempC));
            return (vmbtps);
        }

        //Equation 2.8 
        public  double Vbtps(double vmbtps, double cbtps)
        {
            return (vmbtps / cbtps);
        }

        //Init Test
        public void InitTestConditions(bool useDryAirConditions, uint sensorId, bool sensorCalibrated, double sensorCorrFactor, DateTime calibrationTime,
        bool useInspOrExpirForStart)
        {
            MUseDryAirConditions = useDryAirConditions;
            MSensorId = sensorId;
            MSensorCalibrated = sensorCalibrated;
            MSensorCorrFactor = sensorCorrFactor;
            MCalibrationTime = calibrationTime;

            MFlowVsTime.Clear();
            MVolVsTime.Clear();

            MTestStarted = false;
            MStartTestIndex = 0;

            MSensorZeroed = false;
            MSensorZeroValue = 0;
            MNumSamplesReceived = 0;

            MSensorZeroMinValue = 65536;
            MSensorZeroMaxValue = -65536;

            MUseInspOrExpirForStart = useInspOrExpirForStart;
            MVolumeSum = 0;

            MInitialFlowDetected = false;
            MNumPtsAboveMinFlow = 0;
            MPositiveFlowVolume = 0;
            MStartOfInspiration = 0;
        }

        public void SetSensorZeroed(int tempSample, short sensorZeroValue)
        {
            MSensorZeroed = true;
            MSensorZeroValue = sensorZeroValue;
            SetTemperature(TempInDegC(tempSample));//currently line 607
        }
        public short GetSensorZeroValue()
        {
            return (short)MSensorZeroValue;
        }

        //problems could occur with FlowValue, DeltaVolumeValue, and SensorZeroed
        //SensorZeroed is the only one that can't be referenced as out, should check logic
        //SensorZeroed is suppose to be set by this functoin, but it can't be referenced by out for some reason
        //To fix this, added a statement at line 436 so that sensor zeroed can be set by all loop statements
        //InitTestConditions needs to be called before ProcessData according to original creator

        public void ProcessData(short pressSample, int tempSample, out double flowValue, out double deltaVolumeValue,
        out bool sensorZeroed, bool testStarted)
        {

            //Real code 
            int currentIdx;
            var dFps = MpIirFilter.DoFilter((double)pressSample);
            short filteredPressureSample;

            if (dFps > 0.0)
                dFps = dFps + 0.5;
            else if (dFps < 0.0)
                dFps = dFps - 0.5;

            filteredPressureSample = (short)dFps;

            if (!MSensorZeroed)
            {
                MNumSamplesReceived++;
                if (MNumSamplesReceived > Globals.NumFilterSettlingSamples)
                {
                    MSensorZeroValue += filteredPressureSample;
                    if (filteredPressureSample > MSensorZeroMaxValue)
                    {
                        MSensorZeroMaxValue = filteredPressureSample;
                    }
                    if (filteredPressureSample < MSensorZeroMinValue)
                    {
                        MSensorZeroMinValue = filteredPressureSample;
                    }
                }

                flowValue = 0.0;
                deltaVolumeValue = 0.0;

                if (MNumSamplesReceived == Globals.NumSensorZeroSamples)
                {
                    if (MSensorZeroMaxValue - MSensorZeroMinValue < Globals.SensorZeroMaxPeakToPeak)
                    {
                        MSensorZeroValue = MSensorZeroValue / (Globals.NumSensorZeroSamples - Globals.NumFilterSettlingSamples);
                        MSensorZeroed = true;
                        sensorZeroed = true;
                        SetTemperature(TempInDegC(tempSample));
                        //This line below is used for debugging purposes?(language: c++)
                        //TRACE(TEXT("Sensor Zero Value = %d\n"), MSensorZeroValue);
                    }
                    else
                    {
                        sensorZeroed = false;
                        MSensorZeroMinValue = 65536;
                        MSensorZeroMaxValue = -65536;
                        MNumSamplesReceived = 0;
                        MSensorZeroValue = 0;
                    }
                }
                else
                {
                    sensorZeroed = false;
                    MTestStarted = false;
                }
            }
            else
            {
                sensorZeroed = MSensorZeroed;
                var currentTemp = TempInDegC(tempSample);

                if (MTestStarted == false && testStarted == true)
                {
                    MTestStarted = true;
                    MVolumeSum = 0.0;
                    MStartTestIndex = 0;
                    MTestTemperature.MIsValid = true;
#if DEMO_VERSION
                                        MTestTemperature.MMeasVal = 23.0;
#else
                    MTestTemperature.MMeasVal = currentTemp;
#endif
                    MNumSamplesReceived = 0;
                    flowValue = AdToFlow(filteredPressureSample, (short)MSensorZeroValue, true);
                    InitializeStorage(MVolumeSum, flowValue);
                }

                flowValue = AdToFlow(filteredPressureSample, (short)MSensorZeroValue, true);
                if (MTestStarted)
                {
                    //Set as Grow Line
                    try
                    {
                        MFlowVsTime[MNumSamplesReceived] = flowValue;
                    }
                    catch (ArgumentOutOfRangeException)
                    {
                        Globals.CreateNewList(MFlowVsTime, MNumSamplesReceived, flowValue);
                    }

                    currentIdx = MFlowVsTime.Count();
                    deltaVolumeValue = CalculateDeltaVolume(currentIdx, true);

                    if (MFvcTestStatus != Globals.FvcVolumePlateauReached)
                    {
                        MVolumeSum += deltaVolumeValue;

                        if (flowValue > 0.0)
                        {
                            MPositiveFlowVolume += Math.Abs(deltaVolumeValue);
                        }
                    }

                    try
                    {
                        MVolVsTime[MNumSamplesReceived] = MVolumeSum;
                    }
                    catch (ArgumentOutOfRangeException)
                    {
                        Globals.CreateNewList(MVolVsTime, MNumSamplesReceived, MVolumeSum);
                    }
                    MNumSamplesReceived++;
                    StoreSamples(MVolumeSum, flowValue);

                }
                else
                {
                    deltaVolumeValue = 0.0;
                }
            }
        }

        public double CalculateDeltaVolume(int idx, bool useCorrectionFactor)
        {
            double deltaVol, flowT, flowT1;
            double vmbtps;
            double corrVbtps;
            double sensorCorrFactor;

            var t = idx;

            if (t < 1)
            {
                return 0.0;
            }

            if (useCorrectionFactor)
            {
                sensorCorrFactor = MSensorCorrFactor;
            }
            else
            {
                sensorCorrFactor = 1.0;
            }

            flowT = MFlowVsTime[t];
            flowT1 = MFlowVsTime[t - 1];

            deltaVol = (flowT1 + (flowT - flowT1) / 2) / MSampleRate;

            if (MUseDryAirConditions)
            {
                deltaVol = deltaVol / sensorCorrFactor;
            }
            else
            {
                if (flowT >= 0.0)
                {
                    vmbtps = VMATPS_To_VMBTPS(deltaVol, MTestBarometricPressure.MMeasVal, 33.0);
                    vmbtps = vmbtps * 0.985;
                }
                else
                {
                    vmbtps = VMATPD_To_VMBTPS(deltaVol, MTestBarometricPressure.MMeasVal, MTestTemperature.MMeasVal);
                }
                corrVbtps = Vbtps(vmbtps, sensorCorrFactor);
                deltaVol = corrVbtps;
            }

            return deltaVol;
        }

        //Function needs to be tested
        public void BuildVolFlowFromStoredData(bool useCorrectionFactor)
        {
            int numPoints;
            short deltaPressure;
            var pressValue = MPressStartVal;
            double flowValue, deltaVolumeValue;
            short deltaVol, deltaFlow;
            double volValue;

            if (MVersion <= 1)
            {
                //I hope I am not 1 off.
                numPoints = MPressureData.Count();
                for (var i = 0; i < numPoints; i++)
                {
                    deltaPressure = MPressureData[i];
                    if ((deltaPressure & 0x0080) != 0)
                    {
                        deltaPressure = (short)(deltaPressure | 0xff00);
                    }
                    pressValue += deltaPressure;
                    flowValue = AdToFlow(pressValue, 0, useCorrectionFactor);

                    try
                    {
                        MFlowVsTime[i] = flowValue;
                    }
                    catch (ArgumentOutOfRangeException)
                    {
                        Globals.CreateNewList(MFlowVsTime, i, flowValue);
                    }

                    deltaVolumeValue = CalculateDeltaVolume(i, useCorrectionFactor);
                    MVolumeSum += deltaVolumeValue;


                    try
                    {
                        MVolVsTime[i] = MVolumeSum;
                    }
                    catch (ArgumentOutOfRangeException)
                    {
                        Globals.CreateNewList(MVolVsTime, i, MVolumeSum);
                    }
                }
            }
            else
            {
                numPoints = MVolVsTimeWords.Count();
                volValue = MVolVsTimeStartVal;
                flowValue = MFlowVsTimeStartVal;

                if (MVolVsTime.Count() > 0)
                {
                    MVolVsTime.Clear();
                }
                if (MFlowVsTime.Count() > 0)
                {
                    MFlowVsTime.Clear();
                }

                for (var i = 0; i < numPoints; i++)
                {
                    deltaVol = (short)MVolVsTimeWords[i];
                    volValue = volValue + ((double)deltaVol) / 1000.0;
                    MVolVsTime.Add(volValue);
                    deltaFlow = (short)MFlowVsTimeWords[i];
                    flowValue = flowValue + ((double)deltaFlow) / 1000.0;
                    MFlowVsTime.Add(flowValue);

                }
            }
        }

        public void SetBarometricPressure(double mmHg)
        {
            MTestBarometricPressure.MMeasVal = mmHg;
            MTestBarometricPressure.MIsValid = true;
        }
        public void SetTemperature(double temperature)
        {
            MTestTemperature.MMeasVal = temperature;
            MTestTemperature.MIsValid = true;
        }

        public virtual uint TestTerminationConditions()
        {
            var numPointsSoFar = MFlowVsTime.Count();
            double flowValue;
            int timeSinceFlowDetected;

            if (numPointsSoFar <= 0)
            {
                return 0;
            }

            if (numPointsSoFar > Globals.AbsoluteStoredTestLength)
            {
                return Globals.AbsoluteTimeExceeded;
            }

            if (!MInitialFlowDetected)
            {
                flowValue = MFlowVsTime[numPointsSoFar - 1];
                if (flowValue > Globals.TestStartFlow)
                {
                    MNumPtsAboveMinFlow++;
                    if (MNumPtsAboveMinFlow > Globals.PointsAboveStartFlow)
                    {
                        MInitialFlowDetected = true;
                        MStartTestIndex = numPointsSoFar;
                    }
                }
                else
                {
                    MNumPtsAboveMinFlow = 0;
                }
            }
            else
            {
                timeSinceFlowDetected = numPointsSoFar - MStartTestIndex;
                if (timeSinceFlowDetected > Globals.DefaultNumTestSamples)
                {
                    return Globals.TestTimeExceeded;
                }
            }
            return 0;
        }

        public uint TerminateFVC_VC_CAL(int curIndex, double volThreshold)
        {
            int intervalInSec;
            double pkToPk;
            uint retVal = 0;

            if (curIndex <= 0)
            {
                retVal = 0;
            }

            if (MTestStarted)
            {
                intervalInSec = CalcIntervalInSeconds(curIndex - MStartTestIndex);
                if (intervalInSec >= Globals.MaximumFvcVcCalTest)
                {
                    retVal = Globals.TestTimeExceeded;
                }
                else
                {
                    if (intervalInSec < Globals.TestFileVolumeChangeWindow)
                    {
                        retVal = 0;
                    }
                    else
                    {
                        pkToPk = CalcPeakToPeak(curIndex - Globals.TestFileVolumeChangeWindow * MSampleRate, curIndex, MVolVsTime);
                        if (MVolVsTime[curIndex] == 0)
                        {
                            retVal = Globals.TestTerminatedByOperator;
                        }
                        else if (pkToPk < volThreshold && volThreshold != 0.0)
                        {
                            retVal = Globals.NoVolumeChange;
                        }
                    }
                }
            }
            else
            {
                intervalInSec = CalcIntervalInSeconds(curIndex);
                if (intervalInSec >= Globals.MaximumFvcVcCalTest)
                {
                    retVal = Globals.NoFlowDetected;
                }
            }
            return retVal;
        }

        public uint TerminateMvv(int curIndex)
        {
            int intervalInSec;
            uint retVal = 0;

            if (curIndex <= 0)
                retVal = 0;

            if (MTestStarted)
            {
                intervalInSec = CalcIntervalInSeconds(curIndex - MStartTestIndex);
                if (intervalInSec >= Globals.MaximumMvvTest)
                {
                    retVal = Globals.TestTimeExceeded;
                }
            }
            else
            {
                intervalInSec = CalcIntervalInSeconds(curIndex);
                if (intervalInSec >= Globals.MaximumMvvTest)
                {
                    retVal = Globals.NoFlowDetected;
                }
            }

            return retVal;
        }

        public int CalcIntervalInSeconds(int numSamples)
        {
            return (numSamples / MSampleRate);
        }

        public  double CalcPeakToPeak(int startIdx, int stopIdx, List<double> dataArray)
        {
            var arrayUpper = dataArray.Count();
            double minVal, maxVal, curVal;

            if (stopIdx <= startIdx || arrayUpper <= 0 || stopIdx > arrayUpper)
            {
                return 0.0;
            }

            minVal = float.MaxValue;
            maxVal = -float.MaxValue;

            for (var i = startIdx; i <= stopIdx; i++)
            {
                curVal = dataArray[i];
                if (curVal > maxVal)
                {
                    maxVal = curVal;
                }

                if (curVal < minVal)
                {
                    minVal = curVal;
                }
            }

            return maxVal - minVal;
        }

        public double CalcDifference(int startIdx, int curIdx, List<double> dataArray)
        {
            var arrayUpper = dataArray.Count();
            double startVal, curVal;
            if (curIdx <= startIdx || arrayUpper <= 0 || curIdx > arrayUpper)
            {

            }

            startVal = dataArray[startIdx];
            curVal = dataArray[curIdx];
            return curVal - startVal;
        }


        public  int FindMaxOrMinPeak(int startIdx, int stopIdx, bool findMax, List<double> dataArray)
        {
            var arrayUpper = dataArray.Count();
            int minMaxIdx;
            double testVal, curValue;

            if (startIdx >= stopIdx || stopIdx > arrayUpper || arrayUpper <= 0)
            {
                return (-1);
            }

            if (findMax)
                testVal = -float.MaxValue;
            else
                testVal = float.MaxValue;

            minMaxIdx = -1;
            for (var i = startIdx; i <= stopIdx; i++)
            {
                curValue = dataArray[i];
                if (findMax)
                {
                    if (curValue > testVal)
                    {
                        testVal = curValue;
                        minMaxIdx = i;
                    }
                }
                else
                {
                    if (curValue < testVal)
                    {
                        testVal = curValue;
                        minMaxIdx = i;
                    }
                }
            }

            return minMaxIdx;

        }

        public int RefineMinOrMaxPeak(int maxMinIdx, int numPoints, double epsilon, List<double> dataArray)
        {
            var arrayUpper = dataArray.Count();
            int leftLimitIdx, rightLimitIdx;
            int leftIdx, rightIdx, i;
            double lowerLimit, upperLimit, maxMinVal, curValue;

            leftLimitIdx = maxMinIdx - numPoints;
            if (leftLimitIdx < 0)
                leftLimitIdx = 0;
            rightLimitIdx = maxMinIdx + numPoints;
            if (rightLimitIdx > arrayUpper)
                rightLimitIdx = arrayUpper;

            maxMinVal = dataArray[maxMinIdx];
            lowerLimit = maxMinVal - epsilon;
            upperLimit = maxMinVal + epsilon;

            leftIdx = maxMinIdx;
            curValue = maxMinVal;
            i = leftIdx;
            while (i >= leftLimitIdx && curValue >= lowerLimit && curValue <= upperLimit)
            {
                leftIdx = i;
                curValue = dataArray[i--];
            }

            rightIdx = maxMinIdx;
            i = rightIdx;
            curValue = maxMinVal;

            while (i <= rightLimitIdx && curValue >= lowerLimit && curValue <= upperLimit)
            {
                rightIdx = i;
                curValue = dataArray[i++];
            }
            return ((rightIdx + leftIdx) / 2);
        }

        //I dont take an integer in anymore, it becomes a string for FailTest
        public void FailTest(string idsReason)
        {
            string str;

            MTestAccepted = false;
            MTestAdequate = false;
            MIsBestTest = false;

            //cstring and loadstring
            str = SpirometryApp.ResourceManager.GetString(idsReason);
            if (string.IsNullOrEmpty(MStrReasonNotAdequate))
            {
                MStrReasonNotAdequate = str;
            }
            else
            {
                MStrReasonNotAdequate = MStrReasonNotAdequate + ", " + str;
            }

        }
        //AddReason needs to get a string in as well to replace loadstring
        public void AddReason(string idsReason)
        {
            string strNewReason;

            strNewReason = SpirometryApp.ResourceManager.GetString(idsReason);
            MStrReasonNotAdequate = MStrReasonNotAdequate + ", " + strNewReason;
        }

        //Having made CPredictedValues yet, so assuming that a list is being brought into this function
        public virtual int PercentComplete(List<CPredictedValues> pPredictedValues, bool ispcpMode)
        {
            return 0;
        }

        public void FindBreaths(int startIdx, List<Globals.CBreath> pBreaths, int endIdx)
        {
            var isPositive = true;
            int arrayUpper;
            var pCurrentBreath = new Globals.CBreath();
            var done = false;
            var currentIdx = startIdx;
            var minInterval = Globals.MinBreathInterval;

            if (endIdx != 0)
                arrayUpper = endIdx;   // 28 Oct 00 PHR
            else
                arrayUpper = MFlowVsTime.Count;

            if (MFlowVsTime[currentIdx] > 0)
            {
                pCurrentBreath.MStartIdx = currentIdx;
                pCurrentBreath.MIsFullBreath = false;
            }

            while (!done)
            {
                var nextCrossing = FindZeroCrossing(currentIdx, MFlowVsTime, ref isPositive, arrayUpper);
                if (nextCrossing == -1)
                {
                    done = true;

                    //Need to check if this is equal to pCurrentBreath != NULL when code is refactored
                    if (!pCurrentBreath.BoolCheck())
                    {
                        pCurrentBreath.MStopIdx = MFlowVsTime.Count;
                        //This code never runs through

                        if (pCurrentBreath.MStopIdx - pCurrentBreath.MStopIdx >= minInterval)
                        {
                            pCurrentBreath.MIsFullBreath = false;
                            pCurrentBreath.MVolume = CalculateBreathVolume(pCurrentBreath.MStartIdx, pCurrentBreath.MStopIdx);
                            if (pCurrentBreath.MVolume > Globals.MinimumDeltaVol)
                            {
                                pBreaths.Add(pCurrentBreath);
                            }
                            pCurrentBreath.ReInit();

                        }

                    }

                }
                else
                {
                    if (isPositive)
                    {
                        if (pCurrentBreath.BoolCheck())
                        {
                            pCurrentBreath.MStartIdx = nextCrossing;
                            currentIdx = nextCrossing + minInterval;
                            if (currentIdx > arrayUpper)
                            {
                                currentIdx = arrayUpper;
                            }
                        }
                        else
                        {
                            currentIdx++;
                            if (currentIdx > arrayUpper)
                            {
                                currentIdx = arrayUpper;
                            }
                        }
                    }
                    else
                    {
                        if (!pCurrentBreath.BoolCheck())
                        {
                            pCurrentBreath.MStopIdx = nextCrossing;
                            pCurrentBreath.MIsFullBreath = true;
                            pCurrentBreath.MVolume = CalculateBreathVolume(pCurrentBreath.MStartIdx, pCurrentBreath.MStopIdx);
                            if (pCurrentBreath.MVolume > Globals.MinimumDeltaVol && ((pCurrentBreath.MStopIdx - pCurrentBreath.MStartIdx) > Globals.MinBreathInterval))
                            {
                                pBreaths.Add(pCurrentBreath);
                                pCurrentBreath.ReInit();
                                currentIdx = nextCrossing + minInterval;
                            }
                            else
                            {
                                currentIdx++;
                            }

                            if (currentIdx >= arrayUpper)
                            {
                                done = true;
                            }
                        }
                        else
                        {
                            currentIdx = nextCrossing + minInterval;
                            if (currentIdx >= arrayUpper)
                            {
                                done = true;
                            }
                        }
                    }
                }
            }

            if (!pCurrentBreath.BoolCheck())
            {
                if ((pCurrentBreath.MStopIdx - pCurrentBreath.MStartIdx) > Globals.MinBreathInterval)
                {
                    pCurrentBreath.MVolume = CalculateBreathVolume(pCurrentBreath.MStartIdx, pCurrentBreath.MStopIdx);
                    if (pCurrentBreath.MVolume > Globals.MinimumDeltaVol)
                    {
                        pBreaths.Add(pCurrentBreath);
                        //there are some else statements here comparable to line 1354-1369 and I thought it was redundant, so I deleted it
                    }
                }
                pCurrentBreath.ReInit();
            }

        }

        public  double CalculateBreathVolume(int startIdx, int stopIdx)
        {
            var maxVolIdx = FindMaxOrMinPeak(startIdx, stopIdx, true, MVolVsTime);
            var minVolIdx = FindMaxOrMinPeak(startIdx, stopIdx, false, MVolVsTime);


            if (maxVolIdx >= 0 && maxVolIdx <= MVolVsTime.Count() && minVolIdx >= 0 && minVolIdx <= MVolVsTime.Count())
            {
                var maxVol = MVolVsTime[maxVolIdx];
                var minVol = MVolVsTime[minVolIdx];

                return maxVol - minVol;
            }
            else
            {
                return 0.0;
            }
        }


        public  int FindZeroCrossing(int startIdx, List<double> pArray, ref bool isPositive, int stopIdx)
        {
            var arrayUpper = stopIdx;
            if (startIdx >= arrayUpper || arrayUpper <= 0)
            {
                return -1;
            }

            for (var i = startIdx + 1; i < arrayUpper; i++)
            {
                var prevPt = pArray[i - 1];
                var curPt = pArray[i];

                if (prevPt <= 0.0 && curPt > 0.0)
                {
                    isPositive = true;
                    return i;
                }
                else if (curPt < 0.0 && prevPt >= 0.0)
                {
                    isPositive = false;
                    return i;
                }
            }

            return -1;
        }

        public  void CalcMinAndMax(IEnumerable<double> input, out double minVal, out double maxVal)
        {
            minVal = input.Min();
            maxVal = input.Max();
        }

        public  void CalcMinAndMaxInRange(IEnumerable<double> array, out double minVal, out double maxVal, int startIdx, int stopIdx)
        {
            minVal = array.Skip(startIdx).Take(stopIdx - startIdx).Min();
            maxVal = array.Skip(startIdx).Take(stopIdx  -startIdx).Max();

        }

        public int CalculateTimeScale()
        {
            double testTime = (MVolVsTime.Count + 1.0) / MSampleRate;
            if (testTime > 16.0)
                return 32;
            else if (testTime > 8.0)
                return 16;
            else
                return 8;
        }

        public int CalculateVolumeScale()
        {
            CalcMinAndMax(MVolVsTime,out var volumeMin,out var volumeMax);
            var volumeRange = volumeMax - volumeMin;

            if (volumeRange > 8.0)
                return 16;
            else
                return 8;
        }
        public double TempInDegC(int temperatureSample)
        {
            if (temperatureSample < 0) return 0;

            double retVal = temperatureSample * 2500.0 / 65536;
            retVal = retVal / 20;               // 20 mv per deg. F, now in degrees F
            retVal = (retVal - 32) / 1.8;       // Convert to deg. C

            return retVal;

        }

        public int DegCToTempSampleValue(double tempDegC)
        {
            double temp = tempDegC * 1.8 + 32.0;
            temp = temp * 20;
            temp = temp * 65536 / 2500;
            return (int)temp;
        }

        public  void CalcTimeAndVolumeScales(out int timeScale, out int volumeScale)
        {
            var timeLength = MVolVsTime.Count + 1 - MStartTestIndex;
            timeLength /= MSampleRate;

            if (timeLength > 16.0)
                timeScale = 32;
            else if (timeLength > 8.0)
                timeScale = 16;
            else if (timeLength > 4.0)
                timeScale = 8;
            else
                timeScale = 4;

            CalcMinAndMax(MVolVsTime, out var volMin, out var volMax);
            var volumeRange = volMax - volMin;

            if (volumeRange > 8.0)
                volumeScale = 16;
            else if (volumeRange > 4.0)
                volumeScale = 8;
            else if (volumeRange > 2.0)
                volumeScale = 4;
            else
                volumeScale = 2;

        }

        public void CalcFlowScale(out int FlowScale)
        {
            if (MFlowVsTime.Count <= 0)
            {
                FlowScale = 8;
            }
            else
            {
                CalcMinAndMax(MFlowVsTime, out var flowMin, out var flowMax);
                double flowRange;
                if (Math.Abs(flowMax) > Math.Abs(flowMin))
                    flowRange = Math.Abs(flowMax);
                else
                    flowRange = Math.Abs(flowMin);

                if (flowRange > 8.0)
                    FlowScale = 16;
                else if (flowRange > 4.0)
                    FlowScale = 8;
                else
                    FlowScale = 4;
            }
        }

        public virtual int GetStartFlowVolIdx()
        {
            return 0;
        }

        public virtual int GetEndOfTest()
        {
            return (MVolVsTime.Count());
        }
        //Function called TestForAcceptability suppose to be here but it has no code inside
        //public virtual void TestForAcceptability(ref CSpiroReport pReport, ref CSpiroConfigProtocol pConfigProtocol)
        //{
            
        //}

        //Assume its a list for now
        public virtual double GetPredictedVolume(ref CPredictedValues pPredictedValues)
        {
            return 0.0;
        }

        public virtual double GetMainVolumeMeasurement()
        {
            return 0.0;
        }

        public virtual bool GetDemoModeSample(out short pRawSample)
        {
            var retVal = true;
            MDemoModeSampleCount++;
            pRawSample = 0;
            if (MDemoModeSampleCount > Globals.DefaultNumTestSamples)
                retVal = false;
            return retVal;

        }

        public double FlowToPressure(double flow)
        {
            var upperIdx = MConvTable.NumEntries - 1;
            var done = false;
            var pressure = 0.0;
            var i = 1;
            if (flow < MConvTable.PressFlow[0].FlowVal)
            {
                pressure = MConvTable.PressFlow[0].PressVal;
            }
            else if (flow > MConvTable.PressFlow[upperIdx].FlowVal)
                pressure = MConvTable.PressFlow[upperIdx].PressVal / Globals.InchesPerBit;
            else
            {
                while (!done && i < upperIdx)
                {
                    if (flow >= MConvTable.PressFlow[i].FlowVal && flow <= MConvTable.PressFlow[i + 1].FlowVal)
                    {
                        pressure = MConvTable.PressFlow[i].PressVal + (flow - MConvTable.PressFlow[i].FlowVal) *
                                   (MConvTable.PressFlow[i + 1].PressVal - MConvTable.PressFlow[i].PressVal) /
                                   (MConvTable.PressFlow[i + 1].FlowVal - MConvTable.PressFlow[i].FlowVal);
                        done = true;
                    }
                    i++;
                }
            }

            return pressure;
        }

        public void CalcFlowTimeFromVolTime(List<double> FtArray)
        {
            double SP = 1.0 / MSampleRate;
            int np = (int) (.080 / (2.0 * SP));
            np++;

            double Denom = 0.0;
            double CurrentFlow;
            int LastIdx = MVolVsTime.Count;

            for (int i = 1; i <= np; i++)
            {
                Denom = Denom + 2 * SP * i * i;
            }

            for (int i = 0; i <= LastIdx; i++)
            {
                if (i < np)
                {
                    if (i == 0)
                    {
                        CurrentFlow = 0.0;
                    }
                    else
                    {
                        CurrentFlow = (MVolVsTime[i] - MVolVsTime[i - 1]) / SP;
                    }
                    FtArray.Add(CurrentFlow);
                }
                else if (i > LastIdx - np)
                {
                    FtArray.Add(0.0);
                }
                else
                {
                    CurrentFlow = 0.0;
                    for (int j = -np; j <= np; j++)
                    {
                        CurrentFlow += j * MVolVsTime[i + j];
                    }
                }
            }



        }
    }
}
