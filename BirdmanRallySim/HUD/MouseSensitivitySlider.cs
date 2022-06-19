using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using UnityEngine.UI;


public class MouseSensitivitySlider : MonoBehaviour
{
    private Slider CurrentSlider;

    // Use this for initialization
    void Start()
    {
        CurrentSlider = GetComponent<Slider>();

        if(MyGameManeger.instance.SettingChanged){
            CurrentSlider.value = MyGameManeger.instance.MouseSensitivity*10f;
        }else{
            MyGameManeger.instance.MouseSensitivity = CurrentSlider.value/10f;
        }
    }

    // Update is called once per frame
    void Update()
    {
        
    }

    public void Method()
    {
        MyGameManeger.instance.MouseSensitivity = CurrentSlider.value/10f;
        MyGameManeger.instance.SettingChanged = true;
    }
}