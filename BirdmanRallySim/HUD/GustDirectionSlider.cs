using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using UnityEngine.UI;


public class GustDirectionSlider : MonoBehaviour
{
    private Text scoreText;
    private AerodynamicCalculator script;
    private Slider CurrentSlider;

    // Use this for initialization
    void Start()
    {
        scoreText = GameObject.Find("GustDirection").GetComponent<Text>();
        script = MyGameManeger.instance.Plane.GetComponent<AerodynamicCalculator>();
        CurrentSlider = GetComponent<Slider>();

        if(MyGameManeger.instance.SettingChanged){
            CurrentSlider.value = MyGameManeger.instance.GustDirection/15f;
        }else{
            MyGameManeger.instance.GustDirection = CurrentSlider.value*15f;
        }
        
        scoreText.text = MyGameManeger.instance.GustDirection.ToString("0.000");
    }

    // Update is called once per frame
    void Update()
    {
        
    }

    public void Method()
    {
        MyGameManeger.instance.GustDirection = CurrentSlider.value*15f;
        scoreText.text = MyGameManeger.instance.GustDirection.ToString("0.000");
        MyGameManeger.instance.SettingChanged = true;
    }
}