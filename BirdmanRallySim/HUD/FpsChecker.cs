using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using UnityEngine.UI;
using TMPro;

/// <Summary>
/// シーンのフレームレートを計測して画面に表示するスクリプトです。
/// </Summary>
public class FpsChecker : MonoBehaviour
{
    private Text fpsText;
    private int frameCount;
    private float elapsedTime;

    // Start is called before the first frame update
    void Start()
    {
        fpsText = this.GetComponent<Text>();
    }

    void Update()
    {
        frameCount++;
        elapsedTime += Time.deltaTime;

        if(elapsedTime >= 1.0f)
        {
            float fps = 1.0f*(frameCount/elapsedTime);

            //string fpsRate = "FPS: {fps.ToString("F2")}";
            fpsText.text= "FPS: " + fps.ToString("F2");

            frameCount = 0;
            elapsedTime = 0f;
        }
    }
}