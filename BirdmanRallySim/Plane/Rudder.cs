using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class Rudder : MonoBehaviour
{
    private AerodynamicCalculator script;

    // Start is called before the first frame update
    void Start()
    {
        script = MyGameManeger.instance.Plane.GetComponent<AerodynamicCalculator>();
    }

    // Update is called once per frame
    void FixedUpdate()
    {
        transform.localRotation = Quaternion.AngleAxis(script.dr, Vector3.up);
    }
}
