package miniufo.test.util;

import java.util.regex.Matcher;  
import java.util.regex.Pattern;  


public final class RegUtils{  
	
    /**
     * Email
     */
    public static final String EMAIL = "^\\w+((-\\w+)|(\\.\\w+))*\\@[A-Za-z0-9]+((\\.|-)[A-Za-z0-9]+)*\\.[A-Za-z0-9]+$";
    
    /**
     * Phone number
     * Range of Yidong  : 139 138 137 136 135 134 147 150 151 152 157 158 159 178 182 183 184 187 188
     * Range of Liantong: 130 131 132 155 156 185 186 145 175 176
     * Range of Dianxin : 133 153 177 173 180 181 189
     * Range of virtual : 170 171
     * Reference        : http://www.jihaoba.com/tools/haoduan/
     */
    public static final String PHONE = "^(((13|18)[0-9]{9})|(15[012356789][0-9]{8})|((147|170|171|173|175|176|177|178)[0-9]{8}))$";
    
    /**
     * Chinese characters
     */
    public static final String CHINESE = "^[\\u4E00-\\u9FA5\\uF900-\\uFA2D]+$";
    
    /**
     * Integer
     */
    public static final String INTEGER = "^-?[1-9]\\d*$";
    
    /**
     * Number
     */
    public static final String NUMBER = "^([+-]?)\\d*\\.?\\d+$";
    
    /**
     * Positive number
     */
    public static final String INTEGER_POS = "^[1-9]\\d*$";
    
    /**
     * Floating number
     */
    public static final String FLOAT = "^([+-]?)\\d*\\.\\d+$";
    
    /** 
     * �������� 
     */  
    public static final String FLOAT_POS = "^[1-9]\\d*.\\d*|0.\\d*[1-9]\\d*$";
    /** 
     * �Ƿ�Ϊ���������֣�����0��00��01�����֣� 
     */  
    public static final String INTEGER_WITH_ZERO_POS = "^(([0-9])|([1-9]([0-9]+)))$";
    /** 
     * �Ƿ�Ϊ�������֣���������������������0��00��01�����֣� 
     */  
    public static final String NUMBER_WITH_ZERO = "^((-)?(([0-9])|([1-9]([0-9]+))))$";
    /** 
     * �Ƿ�Ϊ�����ַ��� 
     */  
    public static final String NUMBER_TEXT = "^([0-9]+)$";
    /** 
     * ����(������0�����������������ж��Ƿ��Ҳ�����Ǹ��� 
     */  
    public static final String NUMBER_ALL = "^((-)?(([0-9])|([1-9][0-9]+))(\\.([0-9]+))?)$";
    /** 
     * QQ��5-14λ 
     */  
    public static final String QQ = "^[1-9][0-9]{4,13}$";
    /** 
     * IP��ַ 
     */  
    public static final String IP = "((?:(?:25[0-5]|2[0-4]\\d|[01]?\\d?\\d)\\.){3}(?:25[0-5]|2[0-4]\\d|[01]?\\d?\\d))";
    /** 
     * �ʱ� 
     */  
    public static final String POST_CODE = "[1-9]\\d{5}(?!\\d)";
    /** 
     * ��ͨ���� 
     */  
    public static final String DATE = "^[1-9]\\d{3}-((0[1-9])|(1[0-2]))-((0[1-9])|([1-2][0-9])|(3[0-1]))$";
    /** 
     * �������ڣ������������2�� 
     * ���ڸ�ʽ��2017-10-19 
     * ��2017/10/19 
     * ��2017.10.19 
     * ��2017��10��19�� 
     * ���31����·ݣ�(((01|03|05|07|08|10|12))-((0[1-9])|([1-2][0-9])|(3[0-1]))) 
     * ���30����·ݣ�(((04|06|11))-((0[1-9])|([1-2][0-9])|(30))) 
     * ���29����·ݣ�(02-((0[1-9])|([1-2][0-9]))) 
     */  
    public static final String DATE_COMPLEX = "^(([1-2]\\d{3})(-|/|.|��)((((01|03|05|07|08|10|12))(-|/|.|��)((0[1-9])|([1-2][0-9])|(3[0-1])))|(((04|06|11))(-|/|.|��)((0[1-9])|([1-2][0-9])|(30)))|(02-((0[1-9])|([1-2][0-9]))))(��)?)$";  
      
    /** 
     * ���ӵ����ڣ����������2�� 
     * �������У�������������2�£���ʽ���£�2017-10-19 
     * ������http://www.jb51.net/article/50905.htm�� 
     * ^((?!0000)[0-9]{4}-((0[1-9]|1[0-2])-(0[1-9]|1[0-9]|2[0-8])|(0[13-9]|1[0-2])-(29|30)|(0[13578]|1[02])-31)|([0-9]{2}(0[48]|[2468][048]|[13579][26])|(0[48]|[2468][048]|[13579][26])00)-02-29)$ 
     */  
    public static final String DATE_COMPLEX_LEAP_YEAR = "^(?:(?!0000)[0-9]{4}-(?:(?:0[1-9]|1[0-2])-(?:0[1-9]|1[0-9]|2[0-8])|(?:0[13-9]|1[0-2])-(?:29|30)|(?:0[13578]|1[02])-31)|(?:[0-9]{2}(?:0[48]|[2468][048]|[13579][26])|(?:0[48]|[2468][048]|[13579][26])00)-02-29)$";  
      
    /** 
     * ������ʽУ��,���Ϸ���True 
     * @param regex ������ʽ 
     * @param content У������� 
     * @return 
     * @author lqy 
     */  
    public static boolean isMatch(String regex, CharSequence content){  
        return Pattern.matches(regex, content);  
    }  
      
    /** 
     * У���ֻ����� 
     * @param mobile 
     * @return 
     * @author lqyao 
     */  
    public static final boolean isMoblie(String mobile){  
        boolean flag = false;  
        if (null != mobile && !mobile.trim().equals("") && mobile.trim().length() == 11) {  
            Pattern pattern = Pattern.compile(PHONE);  
            Matcher matcher = pattern.matcher(mobile.trim());  
            flag = matcher.matches();  
        }  
        return flag;  
    }  
      
    /** 
     * У������ 
     * @param value 
     * @return 
     * @author lqyao 
     */  
    public static final boolean isEmail(String value){  
        boolean flag = false;  
        if (null != value && !value.trim().equals("")) {  
            Pattern pattern = Pattern.compile(EMAIL);  
            Matcher matcher = pattern.matcher(value.trim());  
            flag = matcher.matches();  
        }  
        return flag;  
    }  
      
    /** 
     * У������ 
     * @param password 
     * @return ���ȷ��Ϸ���true������Ϊfalse 
     * @author lqyao 
     * @since 2015-09-24 
     */  
    public static final boolean isPassword(String password){  
        boolean flag = false;  
        if (null != password && !password.trim().equals("")) {  
            password = password.trim();  
            if(password.length() >= 6 && password.length() <= 30){  
                return true;  
            }  
        }  
        return flag;  
    }  
      
    /** 
     * У���ֻ���֤�� 
     * @param value 
     * @return ����������ʽ����true�����򷵻�false 
     * @author lqyao 
     * @since 2015-09-24 
     */  
    public static final boolean isPhoneValidateCode(String value){  
        boolean flag = false;  
        if (null != value && !value.trim().equals("")) {  
            Pattern pattern = Pattern.compile("^8\\d{5}$");  
            Matcher matcher = pattern.matcher(value.trim());  
            flag = matcher.matches();  
        }  
        return flag;  
    }  
  
    /** 
     * �ж��Ƿ�ȫ����д��ĸ 
     * @param str 
     * @return 
     */  
    public static boolean isUpperCase(String str){  
        if(str.isEmpty()){  
            return false;  
        }  
        String reg = "^[A-Z]$";  
        return isMatch(reg,str);  
    }  
    /** 
     * �ж��Ƿ�ȫ��Сд��ĸ 
     * @param str 
     * @return 
     */  
    public static boolean isLowercase(String str){  
        if(str.isEmpty()){  
            return false;  
        }  
        String reg = "^[a-z]$";  
        return isMatch(reg,str);  
    }  
      
    public static boolean isIP(String str){  
        if(str.isEmpty()){  
            return false;  
        }  
        return isMatch(IP, str);  
    }  
      
    /** 
     * ���Ϸ���true������30��31��������2�·ݣ����ϸ��У�飩����ʽΪ2017-10-19 
     * @param str 
     * @return 
     */  
    public static boolean isDate(String str){  
        if(str.isEmpty()){  
            return false;  
        }  
        return isMatch(DATE_COMPLEX_LEAP_YEAR, str);  
    }  
      
    /** 
     * ������У�飬����ô�ϸ� 
     * @param str 
     * @return 
     */  
    public static boolean isDateSimple(String str){  
        if(str.isEmpty()){  
            return false;  
        }  
        return isMatch(DATE, str);  
    }  
      
    /** 
     * ����30��31�죬��û�����������2�·� 
     * @param str 
     * @return 
     */  
    public static boolean isDateComplex(String str){  
        if(str.isEmpty()){  
            return false;  
        }  
        return isMatch(DATE_COMPLEX, str);  
    }  
    /** 
     * �ж��Ƿ�Ϊ�����ַ�������0011��10101��01 
     * @param str 
     * @return 
     */  
    public static boolean isNumberText(String str){  
        if(str.isEmpty()){  
            return false;  
        }  
        return isMatch(NUMBER_TEXT, str);  
    }  
    /** 
     * �ж��������͵����֣�����(������0�����������������ж��Ƿ��Ҳ�����Ǹ��� 
     * @param str 
     * @return 
     */  
    public static boolean isNumberAll(String str){  
        if(str.isEmpty()){  
            return false;  
        }  
        return isMatch(NUMBER_ALL, str);  
    }  
      
    /** 
     * �Ƿ�Ϊ���������֣�����0��00��01�����֣� 
     * @param str 
     * @return 
     */  
    public static boolean isIntegerWithZeroPos(String str){  
        if(str.isEmpty()){  
            return false;  
        }  
        return isMatch(INTEGER_WITH_ZERO_POS, str);  
    }  
      
    /** 
     * �Ƿ�Ϊ��������������������������0��00��01�����֣� 
     * @param str 
     * @return 
     */  
    public static boolean isIntegerWithZero(String str){  
        if(str.isEmpty()){  
            return false;  
        }  
        return isMatch(NUMBER_WITH_ZERO, str);  
    }  
      
    /** 
     * ���Ϸ���true,QQ��5-14λ 
     * @param str 
     * @return 
     */  
    public static boolean isQQ(String str){  
        if(str.isEmpty()){  
            return false;  
        }  
        return isMatch(QQ, str);  
    }  
      
      
    public static void main(String[] args) {  
        System.out.println(isMoblie("13430800244"));  
        System.out.println(isMoblie("17730800244"));  
        System.out.println(isMoblie("17630800244"));  
        System.out.println(isMoblie("14730800244"));  
        System.out.println(isMoblie("18330800244"));  
        System.out.println(isMoblie("19330800244"));  
        System.out.println(isMoblie("1333000244"));  
    }
}
