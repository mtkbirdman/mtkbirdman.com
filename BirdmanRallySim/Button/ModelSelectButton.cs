using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using UnityEngine.UI;
using UnityEngine.SceneManagement;

public class ModelSelectButton : MonoBehaviour
{
    private GameObject LoadingText;
    private bool firstPush = false;

    // Start is called before the first frame update
    void Start()
    {
        LoadingText = GameObject.Find("LoadingText");
        LoadingText.SetActive(false);    
    }

    public void OnClick(int number)
    {
        //Debug.Log("Press Start!");
        if (!firstPush)
        {
            switch (number)
            {
                case 0:
                    MyGameManeger.instance.PlaneName = "QX-18";
                    break;
                case 1:
                    MyGameManeger.instance.PlaneName = "QX-19";
                    break;
                case 2:
                    MyGameManeger.instance.PlaneName = "QX-20";
                    break;
                case 3:
                    MyGameManeger.instance.PlaneName = "ARG-2";
                    break;
                case 4:
                    MyGameManeger.instance.PlaneName = "UL01B";
                    break;
                default:
                    break;
            }
            LoadingText.SetActive(true);    
            SceneManager.LoadScene("FlightScene");
        }
    }
}
