using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using UnityEngine.UI;

public class Angle : MonoBehaviour
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
            "\n\n" + script.alpha.ToString("0.000") 
            + "\r\n" + script.beta.ToString("0.000");
    }
}
