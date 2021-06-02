using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using UnityEngine.UI;
using UnityEngine.SceneManagement;

public class Setting : MonoBehaviour
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
            if(Input.GetKey("left shift") || Input.GetKey("right shift")){
                if(Input.GetKeyDown("w")){MyGameManeger.instance.GustMag+=0.5f;RefreshText();};
                if(Input.GetKeyDown("s")){MyGameManeger.instance.GustMag-=0.5f;RefreshText();};
                if(Input.GetKeyDown("a")){MyGameManeger.instance.GustDirection-=15f;RefreshText();};
                if(Input.GetKeyDown("d")){MyGameManeger.instance.GustDirection+=15f;RefreshText();};
            }else{
                if(Input.GetKeyDown("w")){MyGameManeger.instance.Airspeed_TO+=0.1f;RefreshText();};
                if(Input.GetKeyDown("s")){MyGameManeger.instance.Airspeed_TO-=0.1f;RefreshText();};
                if(Input.GetKeyDown("a")){MyGameManeger.instance.alpha_TO+=0.1f;RefreshText();};
                if(Input.GetKeyDown("d")){MyGameManeger.instance.alpha_TO-=0.1f;RefreshText();};
            }
            MyGameManeger.instance.SettingChanged = true;
        }
    }

    void RefreshText()
    {
        scoreText.text = 
            "\n\t" + "SETTING" + "\t\t\t\t"
            + "\n\t"
            + "\n\t" + "AIRSPEED_TO" + "\t\t" + MyGameManeger.instance.Airspeed_TO.ToString("0.000") + "\t\t" + "UP/DOWN -> W/S"
            + "\n\t" + "ALPHA_TO" + "\t\t\t\t" + MyGameManeger.instance.alpha_TO.ToString("0.000") + "\t\t" + "UP/DOWN -> A/D"
            + "\n\t"
            + "\n\t" + "GUST" + "\t\t\t\t\t\t" + MyGameManeger.instance.GustMag.ToString("0.000") + "\t\t" + "UP/DOWN -> Shit + W/S"
            + "\n\t" + "DIRECTION" + "\t\t\t\t" + MyGameManeger.instance.GustDirection.ToString("000") + "\t\t\t" + "LEFT/RIGHT -> Shit + A/D"
            + "\n\t"
            + "\n\t" + "Press R to apply the new value and restart."
            + "\n\t" + "Press E to set the value to default and Restart."
            + "\n\t" + "Press Q to reselect the model.";
    }
}
