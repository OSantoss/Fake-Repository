using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using NUnit.Framework;
using Calculations;

namespace TestCalculations
{
    [TestFixture]
    public class NUnitTest
    {
        [Test]        
        public void TestIirFilter()
        {
            //Arrange
            var filter = new CIirFilter(Globals.SpiroFilt);
            var expected = .0189392;
            //Act
            var actual = filter.DoFilter(Math.Sin(2 * 3.14 * (1 / 15.0)));
            //Assert
            Assert.AreEqual(expected, actual, .0000001);
        }

        [Test]
        public void TestListException()
        {
            //Arrange
            var testList = new List<double>();
            var expected = 20.0;
            var index = 12;
            double actual;
            double actual2;
            //Act
            try
            {
                testList[index] = expected;
            }
            catch(ArgumentOutOfRangeException)
            {
                Globals.CreateNewList(testList, index, expected);
            }

            actual = testList[index];
            
            try
            {
                testList[3] = 10.0;
            }
            catch (ArgumentOutOfRangeException)
            {
                Globals.CreateNewList(testList, 3, 10.0);
            }

            actual2 = testList[3];
            //Assert
            Assert.AreEqual(expected, actual);
            Assert.AreEqual(testList[3], actual2);
        }

        [Test]
        public void CSpiroTestInitializationTest()
        {
            //Arrange
            var CSpiro = new CSpiroTest();
            var expected = DateTime.Now.Day;
            
            //Act
            var actual = CSpiro.MCalibrationTime.Day;
            //Assert
            Assert.AreEqual(expected,actual);

        }

        [Test]
        public void TestInitStorage()
        {
            //Arrange
            var CSpiro = new CSpiroTest();
            //Act
            CSpiro.InitializeStorage(1034,1000);

            var expected = 1034;
            var actual = CSpiro.MPressStartVal;
            //Assert
            Assert.AreEqual(actual,expected);
        }

        [Test]

        public void TestStoreSamples()
        {
            var CSpiro = new CSpiroTest();
            CSpiro.StoreSamples(10.0,10.0);
            var firstVal = CSpiro.MVolVsTimeWords[0];
            var secondVal = CSpiro.MFlowVsTimeWords[0];
            var expected = 10000;

            Assert.AreEqual(firstVal,expected);
        }

        [TestCase(10,10,0.0)]
        //[TestCase(10,20,0.0)]
        public void TestAdToFlow(short firstVal,short secondVal, double expected)
        {
            //Arrange
            var cSpiro = new CSpiroTest();
            //Act
            var actual = cSpiro.AdToFlow(firstVal, secondVal, false);
            //Assert
            Assert.AreEqual(actual,expected);
        }

        [TestCase(10.0, 10.0, 10.0, 8.2111)]
        [TestCase(2.0, 12.0, 10.0, 1.6108)]
        public void TestVMATPS_To_VMBTPS(double vmatpd, double pb, double tempC, double expected)
        {
            //Arrange
            var Cspiro = new CSpiroTest();
            //Act
            var actual = Cspiro.VMATPS_To_VMBTPS(vmatpd, pb, tempC);
            //Assert
            Assert.AreEqual(actual,expected,.0001);
        }

        [TestCase(10.0, 37.73498)]
        [TestCase(20000.0, 37.73498)]
        public void TestPh20(double tempC, double expected)
        {
            //Arrange
            var cSpiro = new CSpiroTest();
            //Act
            var actual = cSpiro.Ph20(tempC);
            //Assert
            Assert.AreEqual(actual,expected,.00001);
        }

        
    }
}
