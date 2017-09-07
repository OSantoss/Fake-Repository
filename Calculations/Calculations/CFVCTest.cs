using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Calculations
{
    public class CFVCTest : CSpiroTest
    {
        public bool MBeeped;
        public int MVersion;
        public Random Factor;
        public double MRandomScaleFactor;
        public CSpireMeas MV_Ext, MFVC, MFEV10, MFEV30, MFEV60, MFEV10_Over_FEV60;
        public CSpireMeas MFEV10_Over_FVC, MFEV30_Over_FVC, m_FEF200_1200, MFEF25, MFEF50, MFEF75, MFEF85;
        public CSpireMeas MFEF25_75, MFEF75_85, MPEF, MPEFLM, MExp_Time, MFIVC, MPIF, MFIV05, MFEV05, MFEV05_Over_FIV05, MFIF50, MFEF50_Over_FIF50;
        public CFVCTest()
        {
            MVersion = Globals.CFVCTest_VERSION;
            Factor = new Random(Environment.TickCount);
            MRandomScaleFactor = Factor.NextDouble();
            MRandomScaleFactor = .98 + MRandomScaleFactor * .02;
            MBeeped = false;

            MV_Ext = new CSpireMeas();
            MFVC = new CSpireMeas();
            MFEV10 = new CSpireMeas();
            MFEV30 = new CSpireMeas();
            MFEV60 = new CSpireMeas();
            MFEV10_Over_FEV60 = new CSpireMeas();
            MFEV10_Over_FVC = new CSpireMeas();
            MFEV30_Over_FVC = new CSpireMeas();
            m_FEF200_1200 =  new CSpireMeas();
            MFEF25 = new CSpireMeas();
            MFEF50 = new CSpireMeas();
            MFEF75 = new CSpireMeas();
            MFEF85 = new CSpireMeas();
            MFEF25_75 = new CSpireMeas();
            MFEF75_85 = new CSpireMeas();
            MPEF = new CSpireMeas();
            MPEFLM = new CSpireMeas();
            MExp_Time = new CSpireMeas();
            MFIVC = new CSpireMeas();
            MPIF = new CSpireMeas();
            MFIV05 = new CSpireMeas();
            MFEV05 = new CSpireMeas();
            MFEV05_Over_FIV05 = new CSpireMeas();
            MFIF50 = new CSpireMeas();
            MFEF50_Over_FIF50 = new CSpireMeas();
        }

        //Serialize Function goes here
        public void Serialize()
        {
            throw new NotImplementedException(
                "Not Implemented because we need to use old c++ code for legacy deserialization");

        }

        public void CalculateMeasurements()
        {
            double MinimumVolume, MinFlow, MaxFlow;

            MTestAccepted = true;
            MTestAdequate = true;

            if (MVolVsTime.Count / MSampleRate < 3)
            {
                FailTest("IDS_INSUFFICIENT_DATA_COLLECTED");
                return;
            }

            var NewStartIndex = MStartTestIndex;
            while (NewStartIndex > 0 && MFlowVsTime[NewStartIndex] > 0.0)
            {
                NewStartIndex--;
            }

            if (NewStartIndex < 0)
            {
                MStartTestIndex = 0;
            }
            else
            {
                MStartTestIndex = NewStartIndex;
            }

            var MaxVolIdx = FindMaxOrMinPeak(MStartTestIndex, MVolVsTime.Count, true, MVolVsTime);
            var MinVolIdx = FindMaxOrMinPeak(0, MStartTestIndex, false, MVolVsTime);

            if (MinVolIdx >= 0)
            {
                MinimumVolume = MVolVsTime[MinVolIdx];
            }
            else
            {
                MinimumVolume = MVolVsTime[0];
                MinVolIdx = 0;
            }

            var MaxIdx = FindMaxOrMinPeak(MinVolIdx, MaxVolIdx, true, MFlowVsTime);
            CalcMinAndMax(MFlowVsTime, out MinFlow, out MaxFlow);
            if (MaxIdx == -1)
            {
                Console.WriteLine("Unable to find peak flow in CFVCTest::CalculateMeasurements");
                return;
            }

            var NumPoints = (int) (MSampleRate * .08);
            var Epsilon = (Globals.FLOW_FULL_SCALE * 2.0 / Globals.AD_COUNTS) * 2.0;
            var RefinedMaxIdx = RefineMinOrMaxPeak(MaxIdx, NumPoints, Epsilon, MFlowVsTime);

            var ArrayUpper = MVolVsTime.Count();
            var LeftIdx = RefinedMaxIdx - NumPoints;
            if (LeftIdx < 0)
            {
                LeftIdx = 0;
            }

            var RightIdx = RefinedMaxIdx + NumPoints;

            if (RightIdx > ArrayUpper)
            {
                RightIdx = ArrayUpper;
            }

            if (RightIdx == LeftIdx)
            {
                Console.WriteLine(
                    "Error: Not enough points to calculate slope for back extrapolation in CFVCTest::CalculateMeasurements");
                return;
            }

            var Slope = (MVolVsTime[RightIdx] - MVolVsTime[LeftIdx]) / (RightIdx - LeftIdx);

            var VolumeAtMaxFlow = MVolVsTime[RefinedMaxIdx] - MinimumVolume;

            if (Slope == 0.0)
            {
                Console.WriteLine("Error: Slope of line for back extrapolation is zero");
                return;
            }

            var Delta = (VolumeAtMaxFlow - MinimumVolume) / Slope;
            var ExpirationStartIdx = RefinedMaxIdx - (int) (Delta + .5);

            var StartIdx = -1 * (VolumeAtMaxFlow - Slope * RefinedMaxIdx) / Slope;
            if (StartIdx < 0.0)
            {
                StartIdx = 0.0;
            }
            ExpirationStartIdx = (int) StartIdx;

            var FractionalPart = StartIdx - ExpirationStartIdx;

            var VolVsTimeUpperBound = MVolVsTime.Count();
            var ExtVol = 0.0;

            if (ExpirationStartIdx < VolVsTimeUpperBound && (ExpirationStartIdx + 1) < VolVsTimeUpperBound)
            {
                ExtVol = MVolVsTime[ExpirationStartIdx] +
                         (MVolVsTime[ExpirationStartIdx + 1] - MVolVsTime[ExpirationStartIdx]) * FractionalPart;
                MV_Ext.MIsValid = true;
                MV_Ext.MMeasVal = ExtVol - MinimumVolume;
            }
            else
            {
                MV_Ext.MIsValid = false;
                MV_Ext.MMeasVal = 0.0;
            }

            if (ExpirationStartIdx < 0)
            {
                Console.WriteLine("Warning: Index of start of inspiration less than 0, correcting");
                ExpirationStartIdx = 0;
            }

            var ExpStopIdx = FindMaxOrMinPeak(MStartTestIndex, MVolVsTime.Count, true, MVolVsTime);

            CalcMinAndMaxInRange(MVolVsTime, out double VolMin, out var VolMax, MStartTestIndex, ExpStopIdx);

            MAbsoluteStartOfExp = MStartTestIndex;

            MFVC.MIsValid = true;
            MFVC.MMeasVal = VolMax - VolMin;

            var Idx_05 = ExpirationStartIdx + MSampleRate / 2;
            if (Idx_05 > ArrayUpper)
            {
                Console.WriteLine("Error: Expiration time less than 1.0 seconds");
                return;
            }

            var Idx_10 = ExpirationStartIdx + MSampleRate;
            if (Idx_10 > ArrayUpper)
            {
                Console.WriteLine("Error: Expiration time less than 1.0 seconds");
                return;
            }

            MFEV10.MIsValid = true;
            MFEV10.MMeasVal = MVolVsTime[Idx_10] - VolMin;

            var Idx_30 = ExpirationStartIdx + MSampleRate * 3;
            if (Idx_30 > ArrayUpper)
            {
                Console.WriteLine("Warning: Expiration time less than 3.0 seconds");
                MFEV30.MIsValid = false;
            }
            else
            {
                MFEV30.MIsValid = true;
                MFEV30.MMeasVal = MVolVsTime[Idx_30] - VolMin;
            }

            var Idx_60 = ExpirationStartIdx + MSampleRate * 6;
            if (Idx_60 <= ArrayUpper)
            {
                MFEV60.MIsValid = true;
                MFEV60.MMeasVal = MVolVsTime[Idx_60] - VolMin;
                if (MFEV60.MMeasVal < 0.0)
                {
                    MFEV60.MIsValid = false; // 26 Sep 00 PHR
                }
                if (MFEV60.MIsValid == true && MFEV60.MMeasVal != 0.0)
                {
                    MFEV10_Over_FEV60.MIsValid = true;
                    MFEV10_Over_FEV60.MMeasVal = (MFEV10.MMeasVal / MFEV60.MMeasVal) * 100.0;
                }
            }

            if (MFVC.MMeasVal != 0.0)
            {
                MFEV10_Over_FVC.MIsValid = true;
                MFEV10_Over_FVC.MMeasVal = (MFEV10.MMeasVal / MFVC.MMeasVal) * 100.0;
                if (MFEV30.MIsValid == true)
                {
                    MFEV30_Over_FVC.MIsValid = true;
                    MFEV30_Over_FVC.MMeasVal = (MFEV30.MMeasVal / MFVC.MMeasVal) * 100.0;
                }
                else MFEV30_Over_FVC.MIsValid = false;
            }

            var Idx_200 = FindVolIdx(ExpirationStartIdx, VolMin + 0.200);
            if (Idx_200 < 0)
            {
                Idx_200 = FindVolIdx(0, VolMin + 0.200);
            }
            var Idx_1200 = FindVolIdx(ExpirationStartIdx, VolMin + 1.2); //lines 2048
            CalcFEF_X_Y(Idx_200, Idx_1200, VolMin + 0.200, VolMin + 1.2, ref m_FEF200_1200, "FEF200-1200");

            var Idx_25Percent =
                CalcVolIdxAndFEF(ref MFEF25, VolMin + MFVC.MMeasVal * 0.25, "FEF25", ExpirationStartIdx);
            var Idx_50Percent =
                CalcVolIdxAndFEF(ref MFEF50, VolMin + MFVC.MMeasVal * 0.50, "FEF50", ExpirationStartIdx);
            var Idx_75Percent =
                CalcVolIdxAndFEF(ref MFEF75, VolMin + MFVC.MMeasVal * 0.75, "FEF75", ExpirationStartIdx);
            var Idx_85Percent =
                CalcVolIdxAndFEF(ref MFEF85, VolMin + MFVC.MMeasVal * 0.85, "FEF85", ExpirationStartIdx);

            CalcFEF_X_Y(Idx_25Percent, Idx_75Percent, MFVC.MMeasVal * 0.25, MFVC.MMeasVal * 0.75, ref MFEF25_75,
                "FEF25-75");
            CalcFEF_X_Y(Idx_75Percent, Idx_85Percent, MFVC.MMeasVal * 0.75, MFVC.MMeasVal * 0.85, ref MFEF75_85,
                "FEF75-85");

            List<double> FtArray = new List<double>();
            CalcFlowTimeFromVolTime(FtArray);
            CalcMinAndMax(FtArray, out MinFlow, out MaxFlow);

            MPEF.MIsValid = true;
            MPEF.MMeasVal = MaxFlow;

            MPEFLM.MIsValid = MPEF.MIsValid;
            MPEFLM.MMeasVal = MPEF.MMeasVal * 60;

            MExp_Time.MIsValid = true;

            if (ExpStopIdx != -1 && ExpStopIdx <= ArrayUpper && ExpStopIdx > ExpirationStartIdx)
            {
                var IdxInterval = (double) (ExpStopIdx - ExpirationStartIdx);
                MExp_Time.MIsValid = true;
                MExp_Time.MMeasVal = IdxInterval / MSampleRate;
            }
            else
            {
                Console.WriteLine(
                    "Warning: Unable to the index of the stop of expiration, skipping inspiration measurements");
                return;
            }

            var RemainingInspiredTime = (double) (MVolVsTime.Count - ExpStopIdx);
            RemainingInspiredTime = RemainingInspiredTime / MSampleRate;

            if (RemainingInspiredTime < 0.5)
            {
                Console.WriteLine("Inspiration period too short. Skipping inspired measurements");
                return;
            }

            CalcMinAndMaxInRange(MFlowVsTime, out MinFlow, out MaxFlow, ExpStopIdx, MFlowVsTime.Count);
            if (MinFlow > -0.5)
            {
                Console.WriteLine("No inspired volume detected. Skipping inspired measurements");
                return;
            }

            var MaxInspiredFlowIdx = FindMaxOrMinPeak(ExpStopIdx, MFlowVsTime.Count,
                false, MFlowVsTime);

            if (MaxInspiredFlowIdx == -1)
            {
                Console.WriteLine("Could not find position of maximum inspired flow. Skipping inspired measurements");
                return;
            }

            MaxInspiredFlowIdx = RefineMinOrMaxPeak(MaxInspiredFlowIdx, NumPoints, Epsilon, MFlowVsTime);
            var StartInspirationIdx = MaxInspiredFlowIdx;

            while (StartInspirationIdx >= ExpStopIdx && MFlowVsTime[StartInspirationIdx] < 0)
            {
                StartInspirationIdx--;
            }

            MStartOfInspiration = StartInspirationIdx;
            MFIVC.MMeasVal = Math.Abs(CalcPeakToPeak(MStartOfInspiration, MVolVsTime.Count, MVolVsTime));

            if (MFIVC.MMeasVal < 0.500)
            {
                MFIVC.MMeasVal = 0.0;
                MFIVC.MIsValid = false;
                return;
            }
            else
            {
                MFIVC.MIsValid = true;
            }

            if (MStartOfInspiration < Idx_60)
            {
                MFEV60.MIsValid = false;
                MFEV10_Over_FEV60.MIsValid = false;
            }

            MPIF.MIsValid = true;
            MPIF.MMeasVal = Math.Abs(MFlowVsTime[MaxInspiredFlowIdx]);

            var FIV_O5Idx = MStartOfInspiration + MSampleRate / 2;

            if (FIV_O5Idx < MVolVsTime.Count)
            {
                MFIV05.MIsValid = true;
                var Fiv05MeasVal = CalcPeakToPeak(MStartOfInspiration, FIV_O5Idx, MVolVsTime);
                MFIV05.MMeasVal = Math.Abs(Fiv05MeasVal);
            }
            else
            {
                MFIV05.MIsValid = false;
            }

            if (MFEV05.MIsValid == true && MFIV05.MMeasVal != 0.0 && MFIV05.MIsValid == true)
            {
                MFEV05_Over_FIV05.MIsValid = true;
                MFEV05_Over_FIV05.MMeasVal = (MFEV05.MMeasVal / MFIV05.MMeasVal) * 100.0;
            }

            var FIF50Idx = FindVolIdx(MStartOfInspiration,MVolVsTime[MStartOfInspiration] - MFIVC.MMeasVal * 0.5);
            if (FIF50Idx != -1)
            {
                MFIF50.MIsValid = true;
                MFIF50.MMeasVal = Math.Abs(MFlowVsTime[FIF50Idx]);
                if (MFIF50.MMeasVal != 0.0 && MFEF50.MIsValid == true)
                {
                    MFEF50_Over_FIF50.MIsValid = true;
                    MFEF50_Over_FIF50.MMeasVal = (MFEF50.MMeasVal / MFIF50.MMeasVal) * 100.0;
                }
            }
        }




        public void CalcFEF_X_Y(int IdxX, int IdxY, double VolX, double VolY, ref CSpireMeas pMeas, string strWarning)
        {
            if (IdxX != -1 && IdxY != -1 && IdxX != IdxY)
            {
                pMeas.MIsValid = true;
                var IdxXInterp = InterpVolIdx(IdxX, VolX);
                var IdxYInterp = InterpVolIdx(IdxY, VolY);
                pMeas.MMeasVal = (VolY - VolX) / ((IdxYInterp - IdxXInterp) / MSampleRate);

            }
            else
            {
                Console.WriteLine("Warning: Index out of range or in error for {0}", strWarning);
                pMeas.MIsValid = false;
                pMeas.MMeasVal = 0.0;
            }
        }

        public double InterpVolIdx(int Idx, double Vol)
        {
            if (!(Vol >= MVolVsTime[Idx] && Vol <= MVolVsTime[Idx + 1]))
                return (Idx);

            var IdxInterp = (double)Idx;
            IdxInterp += (Vol - MVolVsTime[Idx]) / (MVolVsTime[Idx + 1] - MVolVsTime[Idx]);
            return IdxInterp;
        }

        public int CalcVolIdxAndFEF(ref CSpireMeas pMeas, double VolumeValue, string strWarning, int StartIdx)
        {
            var FoundIdx = FindVolIdx(StartIdx, VolumeValue);

            if (FoundIdx != -1)
            {
                pMeas.MIsValid = true;
                pMeas.MMeasVal = MFlowVsTime[FoundIdx];
            }
            else
            {
                Console.WriteLine("Warning: Index out of range for {0}", strWarning);
                pMeas.MIsValid = false;
                pMeas.MMeasVal = 0.0;
            }
            return FoundIdx;
        }

        public int FindVolIdx(int StartIdx, double VolumeValue)
        {
            var ArrayUpper = MVolVsTime.Count;
            if (StartIdx >= ArrayUpper) return -1;

            var PrevVolumeValue = MVolVsTime[StartIdx];
            for (int i = StartIdx; i < ArrayUpper - 1; i++)
            {
                var CurVolumeValue = MVolVsTime[i];
                if (PrevVolumeValue <= VolumeValue && CurVolumeValue >= VolumeValue)
                {
                    if (i == 0) return (0);
                    else return (i - 1);
                }
                else if (PrevVolumeValue >= VolumeValue && CurVolumeValue <= VolumeValue)
                {
                    if (i == 0) return (0);
                    else return (i - 1);
                }
                else
                {
                    PrevVolumeValue = CurVolumeValue;
                }
            }

            return -1;
        }
    }
}
