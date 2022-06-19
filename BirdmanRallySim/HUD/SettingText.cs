using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using UnityEngine.UI;
using UnityEngine.SceneManagement;

public class SettingText : MonoBehaviour
{
    private Text scoreText;
    private AerodynamicCalculator script;
    
    // Start is called before the first frame update
    void Start()
    {
        scoreText = this.GetComponent<Text>();
        script = MyGameManeger.instance.Plane.GetComponent<AerodynamicCalculator>();

        RefreshText();
    }

    // Update is called once per frame
    void Update()
    {
        if(MyGameManeger.instance.SettingActive){
            RefreshText();
        }
    }

    void RefreshText()
    {
        scoreText.text = "\n";
        scoreText.text += MyGameManeger.instance.FlightMode+"\n";

        if(MyGameManeger.instance.CameraSwitch){
            scoreText.text += "FPS"+"\n";
        }else{
            scoreText.text += "TPS"+"\n";
        }
        if(MyGameManeger.instance.HUDActive){
            scoreText.text += "ON"+"\n";
        }else{
            scoreText.text += "OFF"+"\n";
        }
        if(MyGameManeger.instance.MousePitchControl){
            scoreText.text += "Mouse"+"\n";
        }else{
            scoreText.text += "Keyboard"+"\n";
        }
    }
}
