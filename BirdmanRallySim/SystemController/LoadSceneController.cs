using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using UnityEngine.SceneManagement;

public class LoadSceneController : MonoBehaviour
{
    private GameObject Platform;

    // Start is called before the first frame update
    void Start()
    {
        Platform = GameObject.Find("Platform");
        if(MyGameManeger.instance.FlightMode == "TestFlight"){
            Platform.SetActive(false);  
        }
    }

    // Update is called once per frame
    void Update()
    {
        if(Input.GetKeyDown("r") || Input.GetMouseButton(2)){
            Time.timeScale=1f;
            SceneManager.LoadScene("FlightScene");
        }
        /*if(Input.GetKeyDown("s")){
            Time.timeScale=1f;
            SceneManager.LoadScene("ModelSelectScene");
        }*/
        if(Input.GetKeyDown("m")){
            if(MyGameManeger.instance.FlightMode == "BirdmanRally"){
                MyGameManeger.instance.FlightMode = "TestFlight";
            }else if(MyGameManeger.instance.FlightMode == "TestFlight"){
                MyGameManeger.instance.FlightMode = "BirdmanRally";
            }

        }
    }

}
