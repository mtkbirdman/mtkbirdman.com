using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using UnityEngine.UI;


public class TakeoffVelocitySlider : MonoBehaviour
{
    private Text scoreText;
    private AerodynamicCalculator script;
    private Slider CurrentSlider;

    // Use this for initialization
    void Start()
    {
        scoreText = GameObject.Find("TakeoffVelocity").GetComponent<Text>();
        script = MyGameManeger.instance.Plane.GetComponent<AerodynamicCalculator>();
        CurrentSlider = GetComponent<Slider>();

        if(MyGameManeger.instance.SettingChanged){
            CurrentSlider.value = MyGameManeger.instance.Airspeed_TO*10f;
        }else{
            MyGameManeger.instance.Airspeed_TO = CurrentSlider.value*0.1f;
        }
        
        scoreText.text = MyGameManeger.instance.Airspeed_TO.ToString("0.000");
    }

    // Update is called once per frame
    void Update()
    {

    }

    public void Method()
    {
        MyGameManeger.instance.Airspeed_TO = CurrentSlider.value*0.1f;
        scoreText.text = MyGameManeger.instance.Airspeed_TO.ToString("0.000");
        MyGameManeger.instance.SettingChanged = true;
    }
}