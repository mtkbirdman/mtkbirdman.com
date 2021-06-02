using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using UnityEngine.UI;
using UnityEngine.SceneManagement;

public class LoadSceneButton : MonoBehaviour
{
    private bool firstPush = false;

    public void OnClick(string SceneName)
    {
        if (!firstPush)
        {
            SceneManager.LoadScene(SceneName);
            Time.timeScale = 1f;
            firstPush = true;
        }
    }
}
