using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using UnityEngine.SceneManagement;

public class LoadSceanController : MonoBehaviour
{

    // Update is called once per frame
    void Update()
    {
        if(Input.GetKeyDown("r")){
            Time.timeScale=1f;
            SceneManager.LoadScene("FlightScene");
        }
        if(Input.GetKeyDown("e")){
            Time.timeScale=1f;
            MyGameManeger.instance.SettingChanged = false;
            SceneManager.LoadScene("FlightScene");
        }
        if(Input.GetKeyDown("q")){
            Time.timeScale=1f;
            SceneManager.LoadScene("ModelSelectScene");
        }
    }
}
