using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using UnityEngine.UI;

public class Gust : MonoBehaviour
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
    void Update()
    {
        scoreText.text = 
            "\n" + script.LocalGustMag.ToString("0.000") 
            + "\n" + (script.LocalGustDirection).ToString("0");
    }
}
