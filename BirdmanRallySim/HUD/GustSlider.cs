using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using UnityEngine.UI;


public class GustSlider : MonoBehaviour
{
    private Text scoreText;
    private AerodynamicCalculator script;
    private Slider CurrentSlider;

    // Use this for initialization
    void Start()
    {
        scoreText = GameObject.Find("Gust").GetComponent<Text>();
        script = MyGameManeger.instance.Plane.GetComponent<AerodynamicCalculator>();
        CurrentSlider = GetComponent<Slider>();

        if(MyGameManeger.instance.SettingChanged){
            CurrentSlider.value = MyGameManeger.instance.GustMag*10f;
        }else{
            MyGameManeger.instance.GustMag = CurrentSlider.value*0.1f;
        }
        
        scoreText.text = MyGameManeger.instance.GustMag.ToString("0.000");
    }

    // Update is called once per frame
    void Update()
    {
        
    }

    public void Method()
    {
        MyGameManeger.instance.GustMag = CurrentSlider.value*0.1f;
        scoreText.text = MyGameManeger.instance.GustMag.ToString("0.000");
        MyGameManeger.instance.SettingChanged = true;
    }
}