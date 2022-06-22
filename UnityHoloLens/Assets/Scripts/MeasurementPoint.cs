using System;
using UnityEngine;

public class MeasurementPoint
{

    public TransformData transformData { get; set; }
    public double MeasurementData { get; set; }
    public double FrequencyOfMax { get; set; }
    public int MeasurementIndex { get; set; }
    public GameObject Visualization { get; set; }
    //private int timeStamp; add this, (basically the same as MeasurementIndex though)???

    public override string ToString()
    {
        return String.Format("Point {4}, Coordinates: ({0}, {1}, {2}), value: {3}, frequency: {5}", this.transformData.LocalPosition.x, this.transformData.LocalPosition.y, this.transformData.LocalPosition.z, this.MeasurementData, this.MeasurementIndex, this.FrequencyOfMax);
    }

    public override bool Equals(object obj)
    {
        if (obj == null) return false;
        MeasurementPoint objAsMeasurementPoint = obj as MeasurementPoint;
        if (objAsMeasurementPoint == null) return false;
        else return Equals(objAsMeasurementPoint);
    }

    public bool Equals(MeasurementPoint otherPoint)
    {
        if (otherPoint == null) return false;
        return this.MeasurementIndex.Equals(otherPoint.MeasurementIndex);
    }

    public override int GetHashCode()
    {
        return this.MeasurementIndex;
    }
}
