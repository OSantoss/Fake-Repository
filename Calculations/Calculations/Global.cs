using System;
using System.Collections.Generic;
using System.Linq;
using System.Linq.Expressions;
using System.Text;
using System.Threading.Tasks;


namespace Calculations
{
    public class Globals
    {
        //Variables

        public const int CSpiroTestVersion = 3; //Current Version
        public const double DefaultBarometicPressure = 760.0;
        public const int MaxConvEntries = 200;
        public const int CvtMaxEntries = 200;
        public const int SpiroSampleRate = 20;
        public const double InchesPerBit = 0.00030520;
        public const int NumFilterSettlingSamples = 20;
        public const int NumSensorZeroSamples = SpiroSampleRate + NumFilterSettlingSamples;
        public const int SensorZeroMaxPeakToPeak = 10;
        public const uint FvcVolumePlateauReached = 64;
        public const int AbsoluteStoredTestLength = (120 * SpiroSampleRate);
        public const int AbsoluteTimeExceeded = 16;
        public const double TestStartFlow = 0.100;
        public const int PointsAboveStartFlow = 4;
        public const int DefaultNumTestSamples = SpiroSampleRate * 30;
        public const int TestTimeExceeded = 1;
        public const int MaximumFvcVcCalTest = 30;
        public const int TestFileVolumeChangeWindow = 2;
        public const int TestTerminatedByOperator = 8;
        public const int NoVolumeChange = 4;
        public const int NoFlowDetected = 2;
        public const int MaximumMvvTest = 15;
        public const int MinBreathInterval = SpiroSampleRate / 4;
        public const double MinimumDeltaVol = .050;
        public const int NumDemoWaitPoints = 1 * SpiroSampleRate;

        public const int CFVCTest_VERSION = 2;
        public const double FLOW_FULL_SCALE = 14.0;
        public const int AD_COUNTS = 65535;
        public struct DefaultConvType
        {
            public double FlowVal;
            public double PresVal;
        }

        //Default Conversion Table for Calculations
        public static DefaultConvType[] DefConvTable = new DefaultConvType[]
        {
            new DefaultConvType {FlowVal =  -14.0, PresVal = -6.441673 },
            new DefaultConvType {FlowVal = -13.0,  PresVal = -5.608968 },
            new DefaultConvType {FlowVal = -12.0,  PresVal = -4.917260 },
            new DefaultConvType {FlowVal = -11.0,  PresVal = -4.264071 },
            new DefaultConvType {FlowVal = -10.0,  PresVal = -3.632613 },
            new DefaultConvType {FlowVal = -9.0,   PresVal = -3.051023 },
            new DefaultConvType {FlowVal = -8.0,   PresVal = -2.544351 },
            new DefaultConvType {FlowVal = -7.0,   PresVal = -2.064739 },
            new DefaultConvType {FlowVal = -6.0,   PresVal = -1.642342 },
            new DefaultConvType {FlowVal = -5.0,   PresVal = -1.269144 },
            new DefaultConvType {FlowVal = -4.0,   PresVal = -0.875619 },
            new DefaultConvType {FlowVal = -3.0,   PresVal = -0.543806 },
            new DefaultConvType {FlowVal = -2.000, PresVal = -0.276613 },
            new DefaultConvType {FlowVal = -1.500, PresVal = -0.171826 },
            new DefaultConvType {FlowVal = -1.000, PresVal = -0.088508 },
            new DefaultConvType {FlowVal = -0.9,   PresVal = -0.075140 },
            new DefaultConvType {FlowVal = -0.8,   PresVal = -0.062363 },
            new DefaultConvType {FlowVal = -0.7,   PresVal = -0.050968 },
            new DefaultConvType {FlowVal = -0.6,   PresVal = -0.040897 },
            new DefaultConvType {FlowVal = -0.5,   PresVal = -0.031436 },
            new DefaultConvType {FlowVal = -0.4,   PresVal = -0.023297 },
            new DefaultConvType {FlowVal = -0.3,   PresVal = -0.015565 },
            new DefaultConvType {FlowVal = -0.2,   PresVal = -0.009461 },
            new DefaultConvType {FlowVal = -0.1,   PresVal = -0.003968 },
            new DefaultConvType {FlowVal = 0.000,  PresVal =    0.0000 },
            new DefaultConvType {FlowVal = 0.100,  PresVal =  0.004212 },
            new DefaultConvType {FlowVal = 0.200,  PresVal =  0.009583 },
            new DefaultConvType {FlowVal = 0.300,  PresVal =  0.015768 },
            new DefaultConvType {FlowVal = 0.400,  PresVal =  0.023297 },
            new DefaultConvType {FlowVal = 0.500,  PresVal =  0.031130 },
            new DefaultConvType {FlowVal = 0.600,  PresVal =  0.039447 },
            new DefaultConvType {FlowVal = 0.700,  PresVal =  0.049748 },
            new DefaultConvType {FlowVal = 0.800,  PresVal =  0.060023 },
            new DefaultConvType {FlowVal = 0.900,  PresVal =  0.072332 },
            new DefaultConvType {FlowVal = 1.000,  PresVal =  0.084438 },
            new DefaultConvType {FlowVal = 1.500,  PresVal =  0.160332 },
            new DefaultConvType {FlowVal = 2.000,  PresVal =  0.258403 },
            new DefaultConvType {FlowVal = 3.000,  PresVal =  0.516907 },
            new DefaultConvType {FlowVal = 4.000,  PresVal =  0.852087 },
            new DefaultConvType {FlowVal = 5.000,  PresVal =  1.263711 },
            new DefaultConvType {FlowVal = 6.000,  PresVal =  1.715631 },
            new DefaultConvType {FlowVal = 7.000,  PresVal =  2.180451 },
            new DefaultConvType {FlowVal = 8.000,  PresVal =  2.657783 },
            new DefaultConvType {FlowVal = 9.000,  PresVal =  3.188119 },
            new DefaultConvType {FlowVal = 10.00,  PresVal =  3.748874 },
            new DefaultConvType {FlowVal = 11.00,  PresVal =  4.386334 },
            new DefaultConvType {FlowVal = 12.00,  PresVal =  5.091753 },
            new DefaultConvType {FlowVal = 13.00,  PresVal =  5.848039 },
            new DefaultConvType {FlowVal = 14.00,  PresVal =  6.688763 },
        };
        public struct Filter
        {
            public uint Length;       // Size of filter
            public double[] History;           // Pointer to history in filter
            public double[] Coef;              // Pointer to coefficients of filter
        };

        //PRESS_FLOW_TYPE is the exact same struct as DEFAULT_CONV_TYPE
        public struct PressFlowType
        {
            public double PressVal;
            public double FlowVal;
            
        };
        public struct ConvTableType
        {
            public ConvTableType(int maxConvEntries) : this()
            {
                PressFlow = new PressFlowType[CvtMaxEntries];
                NumEntries = maxConvEntries;
            }


            public int NumEntries;
            public PressFlowType[] PressFlow;
            //public PRESS_FLOW_TYPE[]  PressFlow = new PRESS_FLOW_TYPE[200];
        };

        public static double[] SpiroFiltCoef =
        {
            0.046586024440484687,
            -0.328955615800949670,
            0.064522734312176672,
            2.000000000000000000,
            1.000000000000000000,
            -0.453177483675354880,
            0.466513143045749090,
            2.000000000000000000,
            1.000000000000000000
        };

        public static Filter SpiroFilt = new Filter { Length = 2, History = null, Coef = SpiroFiltCoef };

        public  int NumDefaultEntries = 200;

        public class CBreath
        {

            public int MStartIdx;
            public int MStopIdx;
            public bool MIsFullBreath;
            public double MVolume;
            
            public CBreath()
            {
                MStartIdx = 0;
                MStopIdx = 0;
                MIsFullBreath = true;
                MVolume = 0.0;

            }

            public void ReInit()
            {
                MStartIdx = 0;
                MStopIdx = 0;
                MIsFullBreath = true;
                MVolume = 0.0;
            }
            //True is that the class is empty, false is it has different values than original
            public bool BoolCheck()
            {
                if(MStartIdx == 0 && MStopIdx == 0 && MIsFullBreath && MVolume == 0.0)
                {
                    return true;
                }
                else
                {
                    return false;
                }
            }
        }


        //Functions
        //Equates to NUMDEFAULT_ENTRIES
        public static int DefaultEntryLength()
        {
            return DefConvTable.Length;
        }

        public static void CreateNewList(List<double> mainList, int indexOutOfRange, double valueInsert)
        {
            var x = mainList.Count();
            for (var i=x;i<=indexOutOfRange;i++)
            {
                mainList.Add(0);
            }

            mainList[indexOutOfRange] = valueInsert;
        }
    }
}
