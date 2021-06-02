using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using UnityEngine.SceneManagement;

public class Landing : MonoBehaviour
{
    public GameObject button;
    
    // Start is called before the first frame update
    void Start()
    {
        button.SetActive(false);
        MyGameManeger.instance.Landing = false;
    }

    // Update is called once per frame
    void OnCollisionEnter(Collision collision)
    {
        MyGameManeger.instance.Landing = true;
        Time.timeScale = 0f;
        button.SetActive(true);
    }
}
