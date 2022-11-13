using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using UnityEngine.UI;

public class LoadFactor : MonoBehaviour
{
    private Text scoreText;
    private AerodynamicCalculator script;
    
    // Start is called before the first frame update
    void Start()
    {
        scoreText = this.GetComponent<Text>();
        script = MyGameManeger.instance.Plane.GetComponent<AerodynamicCalculator>();
    }

    // Update is called once per frame
    void FixedUpdate()
    {
        scoreText.text = script.Airspeed.ToString("0.000");
    }
}
